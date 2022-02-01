#!/usr/bin/env python3

import click
import pandas as pd



def load_data(novoa_cis, novoa_mms, tombo_mms, tombo_cis, mao_mms, sancar_cis):
    nc = pd.read_csv(novoa_cis, sep='\t').sort_values('CONTEXT')
    nm = pd.read_csv(novoa_mms, sep='\t').sort_values('CONTEXT')
    
    tm = pd.read_csv(tombo_mms, sep='\t').sort_values('CONTEXT')
    tc = pd.read_csv(tombo_cis, sep='\t').sort_values('CONTEXT')

    mm = pd.read_csv(mao_mms, sep='\t')
    aa = list(set(mm['CONTEXT']) ^ set(tm['CONTEXT']))
    bb = pd.DataFrame([aa, [0]*len(aa), [0]*len(aa)]).T
    bb.set_axis(list(mm.columns), axis=1, inplace=True)
    mm = pd.concat([mm, bb], axis=0).reset_index(drop=True).sort_values('CONTEXT')

    sc = pd.read_csv(sancar_cis, sep='\t').sort_values('CONTEXT')

    return nc, nm, tm, tc, mm, sc


@click.command(short_help='Get comparison all and heatmap')
@click.option('-nm', '--novoa_mms', required=True)
@click.option('-nc', '--novoa_cis', required=True)
@click.option('-mm', '--mao_mms', required=True)
@click.option('-sc', '--sancar_cis', required=True)
@click.option('-tm', '--tombo_mms', required=True)
@click.option('-tc', '--tombo_cis', required=True)
@click.option('-o', '--output', required=True)
def main(novoa_cis, novoa_mms, mao_mms, sancar_cis, tombo_mms, 
    tombo_cis, output):

    nc, nm, tm, tc, mm, sc = load_data(
        novoa_cis, novoa_mms, tombo_mms, tombo_cis, mao_mms, sancar_cis
    )




    

if __name__ == '__main__':
    main()