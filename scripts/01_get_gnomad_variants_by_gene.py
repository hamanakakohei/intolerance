#!/usr/bin/env python3

# see: https://gnomad.broadinstitute.org/data#api

import asyncio
from gql import gql, Client
from gql.transport.aiohttp import AIOHTTPTransport
import pandas as pd
import argparse


parser = argparse.ArgumentParser(description="Query gnomAD GraphQL API for variants")
parser.add_argument("--gene_symbol", 		default="RNU4-2", 	help="Gene symbol to query")
parser.add_argument("--reference_genome", 	default="GRCh38", 	help="")
parser.add_argument("--dataset", 		default="gnomad_r4", 	help="")
parser.add_argument("--out_vcf", 		default="output.vcf",   help="")
args = parser.parse_args()


query = gql(f"""
query VariantsInGene {{
  gene(gene_symbol: "{args.gene_symbol}", reference_genome: {args.reference_genome}) {{
    variants(dataset: {args.dataset}) {{
      variant_id
      pos
      genome {{
        ac
        an
        af
      }}
    }}
  }}
}}
""")


async def main():
    transport = AIOHTTPTransport(url="https://gnomad.broadinstitute.org/api")
    client = Client(transport=transport, fetch_schema_from_transport=True)
    result = await client.execute_async(query)
    return(result)


if __name__ == "__main__":
    result = asyncio.run(main())

    data = []
    for v in result['gene']['variants']:
        if v['genome'] is None:
            continue 
        
        data.append({
            'variant_id': v['variant_id'],
            'pos': v['pos'],
            'ac': v['genome']['ac'],
            'an': v['genome']['an'],
            'af': v['genome']['af'],
        })
    
    # データフレームに変換
    # variant_id を分解
    # INFOフィールドを作成（例：AC=xx;AN=yy;AF=zz）
    df = pd.DataFrame(data)
    df[['CHROM', 'POS', 'REF', 'ALT']] = df['variant_id'].str.split('-', expand=True)
    df['CHROM'] = 'chr' + df['CHROM'].astype(str)
    df['POS'] = df['POS'].astype(int)
    df['INFO'] = 'AC=' + df['ac'].astype(str) + ';AN=' + df['an'].astype(str) + ';AF=' + df['af'].astype(str)
    
    # VCF形式のカラムに合わせる
    vcf_df = df[['CHROM', 'POS', 'variant_id', 'REF', 'ALT']]
    vcf_df['QUAL'] = '.'
    vcf_df['FILTER'] = 'PASS'
    vcf_df['INFO'] = df['INFO']
    
    # ヘッダー行を追加
    vcf_header = [
        '##fileformat=VCFv4.2',
        '##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count">',
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    ]
    
    # ファイルに出力
    with open(args.out_vcf, 'w') as f:
        for line in vcf_header:
            f.write(line + '\n')
        vcf_df.to_csv(f, sep='\t', index=False, header=False)

