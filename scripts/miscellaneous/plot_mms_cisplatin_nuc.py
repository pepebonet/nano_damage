
import click
import pandas as pd




# ------------------------------------------------------------------------------
# CLICK
# ------------------------------------------------------------------------------

@click.command(short_help='script to plot damage differences '
'between cisplatin and mms')
@click.option(
    '-dm', '--damage_mms', default='', 
    help='data containing the damage from the mms experiment'
)
@click.option(
    '-dc', '--damage_cisplatin', default='', 
    help='data containing the damage from the cisplatin experiment '
)
@click.option(
    '-o', '--output', default='', help='output folder'
)
def main(damage_mms, damage_cisplatin, output):
    #Obtain data
    mms = pd.read_csv(damage_mms, sep='\t')
    cisplatin = pd.read_csv(damage_cisplatin, sep='\t')


    import pdb;pdb.set_trace()


if __name__ == '__main__':
    main()