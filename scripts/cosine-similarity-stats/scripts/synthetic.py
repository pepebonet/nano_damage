import numpy as np
import pickle
import gzip
import pandas as pd
import os
from sklearn import metrics


class Signatures:

    def __init__(self, signatures_path):

        with gzip.open(signatures_path, 'rb') as f:
            self.table = pickle.load(f)
        self.channel_index = np.array(self.table.index)
        self.etiologies = {'capecitabine': 'SBS17b', 'platinum': 'SBS31', 'flat': 'SBS3',
                           'benzo': 'SBS4', 'poleta': 'SBS9', 'apobec': 'SBS13'}


def signature_recovery(etiology, n_treated, replicate, signatures_path, deconstruction_folder):

    signatures = Signatures(signatures_path)
    case = f'{etiology}.{n_treated}.{replicate}.catalogue.tsv'
    df = pd.read_csv(os.path.join(deconstruction_folder, f'{case}/{case}.processes.tsv'), sep='\t')
    cosines = []
    for col in df.columns:
        reference = signatures.table[signatures.etiologies[etiology]].values
        profile = df[col].values
        cosine = metrics.pairwise.cosine_similarity([reference], [profile])[0][0]
        cosines.append(cosine)
    cosine = max(cosines)
    index = np.argmax(cosines)
    return index, cosine





