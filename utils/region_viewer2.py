import pysam
import gffutils
import argparse
import os
import json
from Bio.Seq import Seq


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate an interactive HTML viewer for pangenome gene models, variants, and AA alignments.")
    parser.add_argument("-v", "--vcf", required=True, help="Input gzipped VCF (tabix-indexed)")
    parser.add_argument("-f", "--fasta", required=True, help="Reference FASTA (indexed)")
    parser.add_argument("-g", "--gff", required=True, help="GFF annotation file")
    parser.add_argument("-c", "--seqname", required=True, help="Chromosome or contig name")
    parser.add_argument("-s", "--start", type=int, required=True, help="Start position (1-based)")
    parser.add_argument("-e", "--stop", type=int, required=True, help="Stop position (1-based)")
    parser.add_argument("-o", "--output", default="pangenome_viewer.html", help="Output HTML file")
    return parser.parse_args()


def init_gff_db(gff_path):
    db_path = gff_path + ".db"
    if not os.path.exists(db_path):
        gffutils.create_db(gff_path, db_path, force=True, keep_order=True, merge_strategy='merge',
                           sort_attribute_values=True)
    return gffutils.FeatureDB(db_path, keep_order=True)


def truncate_seq(seq, max_len=100):
    if len(seq) > max_len:
        return f"{seq[:3]}...{seq[-3:]} ({len(seq)}bp)"
    return seq


def merge_intervals(intervals):
    if not intervals:
        return []
    sorted_inv = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_inv[0]]
    for current in sorted_inv[1:]:
        prev_start, prev_end = merged[-1]
        if current[0] <= prev_end + 1:
            merged[-1] = (prev_start, max(prev_end, current[1]))
        else:
            merged.append(current)
    return merged


def build_sample_aa(fasta_in, seqname, cds_intervals, strand, phase0, sample_vars):
    """Reconstructs sample CDS sequence with VCF variants applied and translates to AA."""
    if not cds_intervals:
        return ""

    min_start = min(c[0] for c in cds_intervals)
    max_end = max(c[1] for c in cds_intervals)

    # Fetch reference span
    ref_span = fasta_in.fetch(seqname, min_start - 1, max_end)

    # Filter variants within CDS span
    span_vars = [v for v in sample_vars if min_start <= v['pos'] <= max_end]
    span_vars.sort(key=lambda x: x['pos'])

    ref_to_sample_idx = {}
    sample_span_chunks = []
    curr_ref_pos = min_start
    sample_pos = 0

    for v in span_vars:
        v_pos = v['pos']
        v_ref = v['ref']
        v_alt = v['alt']
        v_ref_len = len(v_ref)

        if v_pos < curr_ref_pos:
            continue  # Skip overlapping variants

        if v_pos > curr_ref_pos:
            chunk_len = v_pos - curr_ref_pos
            chunk = ref_span[curr_ref_pos - min_start: v_pos - min_start]
            for i in range(chunk_len):
                ref_to_sample_idx[curr_ref_pos + i] = sample_pos + i
            sample_span_chunks.append(chunk)
            sample_pos += chunk_len
            curr_ref_pos = v_pos

        ref_to_sample_idx[v_pos] = sample_pos
        sample_span_chunks.append(v_alt)

        for i in range(1, v_ref_len):
            ref_to_sample_idx[v_pos + i] = sample_pos + len(v_alt)

        sample_pos += len(v_alt)
        curr_ref_pos = v_pos + v_ref_len

    if curr_ref_pos <= max_end:
        chunk_len = max_end - curr_ref_pos + 1
        chunk = ref_span[curr_ref_pos - min_start: max_end - min_start + 1]
        for i in range(chunk_len):
            ref_to_sample_idx[curr_ref_pos + i] = sample_pos + i
        sample_span_chunks.append(chunk)
        sample_pos += chunk_len

    sample_span_str = "".join(sample_span_chunks)

    # Slice CDS intervals
    sorted_cds = sorted(cds_intervals, key=lambda x: x[0])
    sample_cds_chunks = []

    for c_start, c_end in sorted_cds:
        s_idx = ref_to_sample_idx.get(c_start, c_start - min_start)
        if (c_end + 1) in ref_to_sample_idx:
            e_idx = ref_to_sample_idx[c_end + 1]
        else:
            e_idx = ref_to_sample_idx.get(c_end, c_end - min_start) + 1

        chunk = sample_span_str[s_idx:e_idx]
        sample_cds_chunks.append(chunk)

    sample_cds_genomic = "".join(sample_cds_chunks)

    if strand == '-':
        sample_cds_dna = str(Seq(sample_cds_genomic).reverse_complement())
    else:
        sample_cds_dna = sample_cds_genomic

    # Trim initial reading frame phase
    if phase0 and phase0 < len(sample_cds_dna):
        sample_cds_dna = sample_cds_dna[phase0:]

    return str(Seq(sample_cds_dna).translate(to_stop=False))


def extract_gene_models(db, fasta_in, seqname, start, stop, sample_variants, samples):
    genes = []
    for gene in db.region(seqid=seqname, start=start, end=stop, featuretype='gene'):
        raw_exons = [(e.start, e.end) for e in db.children(gene, featuretype='exon')]
        canonical_exons = merge_intervals(raw_exons)

        if gene.strand == '-':
            canonical_exons.sort(key=lambda x: x[0], reverse=True)
        else:
            canonical_exons.sort(key=lambda x: x[0])

        cds_features = list(db.children(gene, featuretype='CDS'))

        # Determine 5' to 3' CDS order to extract phase (column 8)
        cds_info = []
        for c in cds_features:
            try:
                ph = int(c.frame)
            except (ValueError, TypeError):
                ph = 0
            cds_info.append((c.start, c.end, ph))

        if cds_info:
            if gene.strand == '-':
                cds_info_5to3 = sorted(cds_info, key=lambda x: x[0], reverse=True)
            else:
                cds_info_5to3 = sorted(cds_info, key=lambda x: x[0])
            first_phase = cds_info_5to3[0][2]
        else:
            first_phase = 0

        cds_intervals = merge_intervals([(c[0], c[1]) for c in cds_info])

        # Fetch Reference AA
        if cds_intervals:
            sorted_cds = sorted(cds_intervals, key=lambda x: x[0])
            seq_chunks = [fasta_in.fetch(seqname, c[0] - 1, c[1]) for c in sorted_cds]
            genomic_cds = "".join(seq_chunks)

            if gene.strand == '-':
                ref_cds_dna = str(Seq(genomic_cds).reverse_complement())
            else:
                ref_cds_dna = genomic_cds

            if first_phase and first_phase < len(ref_cds_dna):
                ref_cds_dna = ref_cds_dna[first_phase:]

            ref_aa = str(Seq(ref_cds_dna).translate(to_stop=False))
        else:
            ref_aa = ""

        # Fetch Sample AAs
        sample_aa_dict = {}
        if ref_aa:
            for s in samples:
                sample_aa_dict[s] = build_sample_aa(fasta_in, seqname, cds_intervals, gene.strand, first_phase,
                                                    sample_variants[s])

        gene_data = {
            "id": gene.id or gene.attributes.get('Name', ['Unknown'])[0],
            "start": gene.start,
            "end": gene.end,
            "strand": gene.strand,
            "exons": [{"start": e[0], "end": e[1]} for e in canonical_exons],
            "cds": [{"start": c[0], "end": c[1]} for c in cds_intervals],
            "ref_aa": ref_aa,
            "sample_aa": sample_aa_dict
        }
        genes.append(gene_data)
    return genes


def extract_variants(vcf_path, seqname, start, stop):
    vcf_in = pysam.VariantFile(vcf_path)
    samples = list(vcf_in.header.samples)
    sample_variants = {sample: [] for sample in samples}

    try:
        records = list(vcf_in.fetch(seqname, start, stop))
    except ValueError:
        return sample_variants, samples

    for rec in records:
        var_type = "SNP" if len(rec.ref) == 1 and all(len(a) == 1 for a in rec.alts) else "INDEL/SV"
        v_start = rec.pos

        for sample in samples:
            gt = rec.samples[sample]['GT']
            if gt is not None and any(g > 0 for g in gt if g is not None):
                alt_idx = [g - 1 for g in gt if g is not None and g > 0][0]
                alt_allele = rec.alts[alt_idx]

                sample_variants[sample].append({
                    "pos": v_start,
                    "type": var_type,
                    "ref": rec.ref,
                    "alt": alt_allele,
                    "ref_display": truncate_seq(rec.ref),
                    "alt_display": truncate_seq(alt_allele),
                    "len": len(rec.ref)
                })

    return sample_variants, samples


def generate_html(seqname, start, stop, genes, sample_variants, samples, output_file):
    html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Pangenome Viewer: {seqname}:{start}-{stop}</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; background-color: #f8f9fa; margin: 0; padding: 20px; display: flex; flex-direction: column; gap: 20px; }}
        .container {{ width: 100%; box-sizing: border-box; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); overflow-x: auto; }}
        h1, h2 {{ font-size: 20px; color: #333; margin-top: 0; }}
        .track-label {{ font-size: 12px; font-weight: bold; fill: #555; }}
        .axis-label {{ font-size: 10px; fill: #888; }}
        .exon-pos {{ fill: #28a745; stroke: #1e7e34; stroke-width: 1; cursor: pointer; }}
        .exon-neg {{ fill: #dc3545; stroke: #a71d2a; stroke-width: 1; cursor: pointer; }}
        .intron {{ stroke: #6c757d; stroke-width: 2; cursor: pointer; }}
        .variant-snp {{ fill: #17a2b8; fill-opacity: 0.8; }}
        .variant-indel {{ fill: #fd7e14; fill-opacity: 0.8; }}
        .tooltip {{ position: absolute; padding: 8px; font-size: 12px; background: #333; color: #fff; border-radius: 4px; pointer-events: none; opacity: 0; transition: opacity 0.2s; z-index: 10; white-space: pre-wrap; }}
        #detail-view {{ display: none; background: #fff; padding: 20px; border-radius: 8px; border: 2px solid #007bff; }}
        .detail-variant-list {{ max-height: 350px; overflow-y: auto; font-size: 12px; }}
        .badge {{ display: inline-block; padding: 3px 8px; font-size: 11px; font-weight: 600; border-radius: 4px; }}
        .badge-exon {{ background-color: #d4edda; color: #155724; border: 1px solid #c3e6cb; }}
        .badge-boundary {{ background-color: #fff3cd; color: #856404; border: 1px solid #ffeeba; }}
        .badge-intron {{ background-color: #e2e3e5; color: #383d41; border: 1px solid #d6d8db; }}
        .filter-box {{ margin: 15px 0; display: flex; gap: 20px; align-items: center; background: #f1f3f5; padding: 10px 15px; border-radius: 6px; font-size: 13px; }}

        /* AA Grid Alignment Styles */
        .aa-grid {{ font-family: monospace; font-size: 13px; white-space: pre; overflow-x: auto; background: #272822; color: #f8f8f2; padding: 15px; border-radius: 6px; margin-top: 15px; }}
        .aa-row {{ display: flex; align-items: center; margin-bottom: 4px; }}
        .aa-header {{ width: 130px; font-weight: bold; color: #a6e22e; flex-shrink: 0; }}
        .aa-seq {{ display: flex; letter-spacing: 2px; }}
        .aa-char {{ width: 14px; text-align: center; display: inline-block; }}
        .aa-mut {{ background-color: #f92672; color: #fff; font-weight: bold; border-radius: 2px; }}
        .aa-syn {{ color: #75715e; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Pangenome Region: {seqname}:{start}-{stop}</h1>
        <div id="tooltip" class="tooltip"></div>
        <svg id="viewer" width="100%" height="400"></svg>
    </div>

    <div class="container" id="detail-view">
        <h2 id="detail-title">Gene Details</h2>
        <div id="detail-content"></div>
    </div>

    <script>
        const start = {start};
        const stop = {stop};
        const length = stop - start;
        const genes = {json.dumps(genes)};
        const sampleVariants = {json.dumps(sample_variants)};
        const samples = {json.dumps(samples)};

        const svg = document.getElementById('viewer');
        const tooltip = document.getElementById('tooltip');
        const detailView = document.getElementById('detail-view');
        const detailTitle = document.getElementById('detail-title');
        const detailContent = document.getElementById('detail-content');

        const svgWidth = 1800;
        svg.setAttribute('viewBox', `0 0 ${{svgWidth}} 400`);
        const margin = {{ left: 120, right: 20, top: 40, bottom: 20 }};
        const drawWidth = svgWidth - margin.left - margin.right;
        const trackHeight = 40;

        function scale(pos) {{
            return margin.left + ((pos - start) / length) * drawWidth;
        }}

        const geneRows = [];
        genes.forEach(gene => {{
            let placed = false;
            for (let i = 0; i < geneRows.length; i++) {{
                let overlap = geneRows[i].some(existing =>
                    (gene.start < existing.end + 1000) && (gene.end > existing.start - 1000)
                );
                if (!overlap) {{
                    geneRows[i].push(gene);
                    gene.row = i;
                    placed = true;
                    break;
                }}
            }}
            if (!placed) {{
                gene.row = geneRows.length;
                geneRows.push([gene]);
            }}
        }});

        let svgContent = '';

        // Axis
        svgContent += `<line x1="${{margin.left}}" y1="20" x2="${{svgWidth - margin.right}}" y2="20" stroke="#ccc" stroke-width="2"/>`;
        for (let i = 0; i <= 5; i++) {{
            const pos = Math.round(start + (length * (i / 5)));
            const x = scale(pos);
            svgContent += `<line x1="${{x}}" y1="15" x2="${{x}}" y2="25" stroke="#ccc" stroke-width="2"/>`;
            svgContent += `<text x="${{x}}" y="10" class="axis-label" text-anchor="middle">${{pos.toLocaleString()}}</text>`;
        }}

        let currentY = margin.top;
        svgContent += `<text x="10" y="${{currentY + 25}}" class="track-label">Gene Models</text>`;

        // Gene Models
        genes.forEach((gene, idx) => {{
            const rowY = currentY + (gene.row * 30);
            const gStart = Math.max(start, gene.start);
            const gEnd = Math.min(stop, gene.end);
            const x1 = scale(gStart);
            const x2 = scale(gEnd);
            const exonClass = gene.strand === '+' ? 'exon-pos' : 'exon-neg';

            svgContent += `<g onclick="showGeneDetails(${{idx}})">`;
            svgContent += `<line x1="${{x1}}" y1="${{rowY + 20}}" x2="${{x2}}" y2="${{rowY + 20}}" class="intron"/>`;

            gene.exons.forEach(exon => {{
                if (exon.end >= start && exon.start <= stop) {{
                    const ex = scale(Math.max(start, exon.start));
                    const ew = Math.max(2, scale(Math.min(stop, exon.end)) - ex);
                    svgContent += `<rect x="${{ex}}" y="${{rowY + 12}}" width="${{ew}}" height="16" class="${{exonClass}}"
                        onmouseover="showTooltip(event, '${{gene.id}} (${{gene.strand}})')" onmouseout="hideTooltip()"/>`;
                }}
            }});
            svgContent += `<text x="${{x1}}" y="${{rowY + 8}}" font-size="10" fill="#333">${{gene.id}}</text>`;
            svgContent += `</g>`;
        }});

        currentY += (geneRows.length * 30) + 20;

        // Sample Tracks
        samples.forEach(sample => {{
            svgContent += `<text x="10" y="${{currentY + 20}}" class="track-label">${{sample}}</text>`;
            svgContent += `<line x1="${{margin.left}}" y1="${{currentY + 15}}" x2="${{svgWidth - margin.right}}" y2="${{currentY + 15}}" stroke="#eee" stroke-width="1"/>`;

            if (sampleVariants[sample]) {{
                sampleVariants[sample].forEach(v => {{
                    if (v.pos >= start && v.pos <= stop) {{
                        const vx = scale(v.pos);
                        const vClass = v.type === 'SNP' ? 'variant-snp' : 'variant-indel';
                        const w = v.type === 'SNP' ? 3 : Math.max(4, scale(v.pos + v.len) - vx);
                        svgContent += `<rect x="${{vx - (w/2)}}" y="${{currentY + 10}}" width="${{w}}" height="10" class="${{vClass}}"
                            onmouseover="showTooltip(event, '${{v.type}}: ${{v.pos}}')" onmouseout="hideTooltip()"/>`;
                    }}
                }});
            }}
            currentY += Math.max(25, trackHeight - 15);
        }});

        svg.innerHTML = svgContent;
        svg.setAttribute('height', currentY + margin.bottom);
        svg.setAttribute('viewBox', `0 0 ${{svgWidth}} ${{currentY + margin.bottom}}`);

        function showTooltip(evt, text) {{
            tooltip.innerText = text;
            tooltip.style.opacity = 1;
            tooltip.style.left = (evt.pageX + 10) + 'px';
            tooltip.style.top = (evt.pageY + 10) + 'px';
        }}

        function hideTooltip() {{ tooltip.style.opacity = 0; }}

        window.showGeneDetails = function(geneIdx) {{
            const gene = genes[geneIdx];
            detailView.style.display = 'block';
            detailTitle.innerText = `Detailed Breakdown: ${{gene.id}} (${{gene.strand}} strand)`;

            function getLocalization(v) {{
                const vStart = v.pos;
                const vEnd = v.pos + (v.len || 1) - 1;

                let overlapsExon = false;
                let overlapsIntron = false;
                let hitExonNum = null;

                gene.exons.forEach((exon, idx) => {{
                    const eNum = idx + 1;
                    const intersects = (vStart <= exon.end && vEnd >= exon.start);
                    if (intersects) {{
                        overlapsExon = true;
                        hitExonNum = eNum;
                        if (vStart < exon.start || vEnd > exon.end) {{
                            overlapsIntron = true;
                        }}
                    }}
                }});

                if (!overlapsExon) {{
                    return {{ label: 'Intron', category: 'intron' }};
                }} else if (overlapsExon && overlapsIntron) {{
                    return {{ label: `Exon-Intron Boundary (Exon ${{hitExonNum}})`, category: 'boundary' }};
                }} else {{
                    return {{ label: `Exon ${{hitExonNum}}`, category: 'exon' }};
                }}
            }}

            let html = `<div><strong>Location:</strong> ${{gene.start.toLocaleString()}} - ${{gene.end.toLocaleString()}}</div>`;

            html += `
            <div class="filter-box">
                <strong>Filter Localization:</strong>
                <label style="cursor:pointer;"><input type="checkbox" id="filter-exon" checked onchange="filterVariantTable()"> Exons</label>
                <label style="cursor:pointer;"><input type="checkbox" id="filter-boundary" checked onchange="filterVariantTable()"> Exon-Intron Boundaries</label>
                <label style="cursor:pointer;"><input type="checkbox" id="filter-intron" checked onchange="filterVariantTable()"> Introns</label>
            </div>
            `;

            html += `<div class="detail-variant-list">`;
            html += `<table style="width:100%; text-align:left; border-collapse: collapse; margin-top: 10px;">`;
            html += `<tr style="border-bottom: 2px solid #dee2e6;">
                <th style="padding: 8px;">Sample</th>
                <th style="padding: 8px;">Type</th>
                <th style="padding: 8px;">Position</th>
                <th style="padding: 8px;">Localization</th>
                <th style="padding: 8px;">Ref &rarr; Alt</th>
            </tr>`;

            let variantCount = 0;
            samples.forEach(sample => {{
                if (sampleVariants[sample]) {{
                    sampleVariants[sample].forEach(v => {{
                        if (v.pos >= gene.start && v.pos <= gene.end) {{
                            variantCount++;
                            const loc = getLocalization(v);
                            html += `<tr class="variant-row" data-category="${{loc.category}}" style="border-bottom: 1px solid #e9ecef;">
                                <td style="padding: 8px;">${{sample}}</td>
                                <td style="padding: 8px;">${{v.type}}</td>
                                <td style="padding: 8px;">${{v.pos.toLocaleString()}}</td>
                                <td style="padding: 8px;"><span class="badge badge-${{loc.category}}">${{loc.label}}</span></td>
                                <td style="padding: 8px; font-family: monospace;">${{v.ref_display}} &rarr; ${{v.alt_display}}</td>
                            </tr>`;
                        }}
                    }});
                }}
            }});

            if (variantCount === 0) {{
                html += `<tr><td colspan="5" style="padding: 12px; text-align: center; color: #6c757d;">No variants found in this gene model.</td></tr>`;
            }}
            html += `</table></div>`;

            // Dynamic Amino Acid Alignment View
            if (gene.ref_aa) {{
                html += `<h3 style="margin-top: 25px; font-size: 16px;">Coding Sequence (CDS) Amino Acid Alignment</h3>`;
                html += `<div class="aa-grid">`;

                // Reference AA Row
                html += `<div class="aa-row"><div class="aa-header">Reference</div><div class="aa-seq">`;
                for (let i = 0; i < gene.ref_aa.length; i++) {{
                    html += `<span class="aa-char">${{gene.ref_aa[i]}}</span>`;
                }}
                html += `</div></div>`;

                // Reconstructed Sample AA Rows
                samples.forEach(sample => {{
                    const sampleAA = gene.sample_aa[sample] || "";
                    const maxLen = Math.max(gene.ref_aa.length, sampleAA.length);

                    html += `<div class="aa-row"><div class="aa-header">${{sample}}</div><div class="aa-seq">`;
                    for (let i = 0; i < maxLen; i++) {{
                        const rChar = gene.ref_aa[i] || "-";
                        const sChar = sampleAA[i] || "-";

                        if (sChar === rChar) {{
                            html += `<span class="aa-char aa-syn">${{sChar}}</span>`;
                        }} else {{
                            html += `<span class="aa-char aa-mut" title="AA Pos ${{i+1}}: Ref ${{rChar}} &rarr; ${{sChar}}">${{sChar}}</span>`;
                        }}
                    }}
                    html += `</div></div>`;
                }});

                html += `</div>`;
            }}

            detailContent.innerHTML = html;
            detailView.scrollIntoView({{ behavior: 'smooth' }});
        }};

        window.filterVariantTable = function() {{
            const showExon = document.getElementById('filter-exon').checked;
            const showBoundary = document.getElementById('filter-boundary').checked;
            const showIntron = document.getElementById('filter-intron').checked;

            const rows = document.querySelectorAll('.variant-row');
            rows.forEach(row => {{
                const cat = row.getAttribute('data-category');
                if ((cat === 'exon' && showExon) ||
                    (cat === 'boundary' && showBoundary) ||
                    (cat === 'intron' && showIntron)) {{
                    row.style.display = '';
                }} else {{
                    row.style.display = 'none';
                }}
            }});
        }};
    </script>
</body>
</html>"""

    with open(output_file, 'w') as f:
        f.write(html_template)
    print(f"Interactive viewer written to {output_file}")


def main():
    args = parse_args()

    print("Extracting pangenome variants...")
    sample_variants, samples = extract_variants(args.vcf, args.seqname, args.start, args.stop)

    print("Parsing GFF annotations & translating CDS regions...")
    fasta_in = pysam.FastaFile(args.fasta)
    db = init_gff_db(args.gff)
    genes = extract_gene_models(db, fasta_in, args.seqname, args.start, args.stop, sample_variants, samples)

    print("Generating HTML viewer...")
    generate_html(args.seqname, args.start, args.stop, genes, sample_variants, samples, args.output)


if __name__ == "__main__":
    main()