import pysam
import gffutils
import argparse
import os
import json


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate an interactive HTML viewer for pangenome gene models and sample variants.")
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


def extract_gene_models(db, seqname, start, stop):
    genes = []
    # Fetch all genes overlapping the window
    for gene in db.region(seqid=seqname, start=start, end=stop, featuretype='gene'):
        gene_data = {
            "id": gene.id or gene.attributes.get('Name', ['Unknown'])[0],
            "start": gene.start,
            "end": gene.end,
            "strand": gene.strand,
            "exons": []
        }
        # Get child exons for the intron/exon structure
        for exon in db.children(gene, featuretype='exon'):
            gene_data["exons"].append({
                "start": exon.start,
                "end": exon.end
            })
        genes.append(gene_data)
    return genes


def extract_variants(vcf_path, seqname, start, stop):
    vcf_in = pysam.VariantFile(vcf_path)
    samples = list(vcf_in.header.samples)

    sample_variants = {sample: [] for sample in samples}

    try:
        records = vcf_in.fetch(seqname, start, stop)
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
        .detail-variant-list {{ max-height: 300px; overflow-y: auto; font-family: monospace; font-size: 12px; }}
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
        <div id="detail-content" class="detail-variant-list"></div>
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

        const svgWidth = 1200; // Fixed width for predictable rendering
        svg.setAttribute('viewBox', `0 0 ${{svgWidth}} 400`);
        const margin = {{ left: 120, right: 20, top: 40, bottom: 20 }};
        const drawWidth = svgWidth - margin.left - margin.right;
        const trackHeight = 40;

        function scale(pos) {{
            return margin.left + ((pos - start) / length) * drawWidth;
        }}

        // Track Packing Algorithm for Genes
        const geneRows = [];
        genes.forEach(gene => {{
            let placed = false;
            for (let i = 0; i < geneRows.length; i++) {{
                let overlap = geneRows[i].some(existing =>
                    (gene.start < existing.end + 1000) && (gene.end > existing.start - 1000) // 1kb padding
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

        // Draw Genes tightly packed
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

        // Draw Sample Variants
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

            let html = `<strong>Location:</strong> ${{gene.start.toLocaleString()}} - ${{gene.end.toLocaleString()}}<br><br>`;
            html += `<strong>Variants Mapping to Gene Region:</strong><br><table style="width:100%; text-align:left;">`;
            html += `<tr><th>Sample</th><th>Type</th><th>Position</th><th>Ref &rarr; Alt</th></tr>`;

            samples.forEach(sample => {{
                if(sampleVariants[sample]) {{
                    sampleVariants[sample].forEach(v => {{
                        if (v.pos >= gene.start && v.pos <= gene.end) {{
                            html += `<tr>
                                <td>${{sample}}</td>
                                <td>${{v.type}}</td>
                                <td>${{v.pos.toLocaleString()}}</td>
                                <td>${{v.ref}} &rarr; ${{v.alt}}</td>
                            </tr>`;
                        }}
                    }});
                }}
            }});

            html += `</table>`;
            detailContent.innerHTML = html;
            detailView.scrollIntoView({{ behavior: 'smooth' }});
        }};
    </script>
</body>
</html>"""

    with open(output_file, 'w') as f:
        f.write(html_template)
    print(f"Interactive viewer written to {output_file}")

def main():
    args = parse_args()

    print("Parsing GFF annotations...")
    db = init_gff_db(args.gff)
    genes = extract_gene_models(db, args.seqname, args.start, args.stop)

    print("Extracting pangenome variants...")
    sample_variants, samples = extract_variants(args.vcf, args.seqname, args.start, args.stop)

    print("Generating HTML viewer...")
    generate_html(args.seqname, args.start, args.stop, genes, sample_variants, samples, args.output)


if __name__ == "__main__":
    main()