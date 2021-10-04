import pandas as pd
import numpy as np
import click

class Design:
    def __init__(self, file, control, pairwise=False):
        self.file = file
        self.sample_sheet = pd.read_csv(file, sep = "\t")
        self.control_name = control
        self.conditions = list(set(self.sample_sheet.condition))
        self.samples_dict = None
        self.replicates = []
        self.design_matrix = None

        if pairwise != False:
            self.conditions = pairwise


    def createSamplesSummary(self):
        self.samples_dict = dict.fromkeys(list(self.conditions))
        for col in self.conditions:
            self.samples_dict[col] = dict.fromkeys(['isControl', 'replicates'])
            self.samples_dict[col]['replicates'] = []
            if col == self.control_name:
                self.samples_dict[col]['isControl'] = 1
            else:
                self.samples_dict[col]['isControl'] = 0
            self.samples_dict[col]['replicates']  = self.sample_sheet.loc[self.sample_sheet.condition == col]['replicate'].values

    def getAllReplicates(self):
        for condition in self.conditions:
            self.replicates.extend(self.samples_dict[condition]['replicates'])



    def fillConditions(self):
        for replicate in self.replicates:
            if replicate in self.samples_dict[self.control_name]['replicates']:
                for col in self.conditions:
                    self.design_matrix.at[replicate, col] = 0


    def fillNan(self):
        for condition in self.conditions:
            for replicate in self.samples_dict[condition]['replicates']:
                if condition != self.control_name:
                    self.design_matrix.at[replicate, condition] = 1
        self.design_matrix = self.design_matrix.replace(np.nan,0)

    def get_empty_design_matrix(self):
        dm_cols_name = ['Samples', 'baseline']
        dm_cols_name.extend(list(self.conditions))
        self.design_matrix = pd.DataFrame(columns=dm_cols_name)
        self.design_matrix['Samples'] = self.replicates
        self.design_matrix['baseline'] = 1
        self.design_matrix.index = self.replicates


    def create_design_matrix(self):
        self.createSamplesSummary()
        self.getAllReplicates()
        self.get_empty_design_matrix()
        self.fillConditions()
        self.fillNan()


def get_all_pairwise_comparaisons(samples_list):
    def compareTuples(t1, t2):
        comparaison = set(t1) & set(t2)
        return (len(comparaison) == len(t1)) & (len(comparaison) == len(t2))


    def tupleInList(t, comparaisons):
        inList = False
        for comparaison in comparaisons:
            if compareTuples(t, comparaison):
                inList = True
                break
        return inList


    def pairwiseComparaisons(samples_list):
        comparaisons = []
        for i in range(0, len(samples_list)):
            for j in range(0, len(samples_list)):
                if (i != j) and not (tupleInList((i, j), comparaisons)):
                    comparaisons.append((i, j))
        return comparaisons

    comparaisons = pairwiseComparaisons(samples_list)
    return comparaisons

@click.command()
@click.option('--file', default=1, help='Samples File.', required=True, type=str)
@click.option('--directory', default=1, help='Output directory.', required=True, type=str)
@click.option('--control', default=1, help='Control name.', required=True, type=str)
@click.option('--treatment', default=1, help='Treatment name.', required=True, type=str)
@click.option('--dryrun', default=False, help='Dry run.', required=True, type=bool)
def readSamples(file, directory, control, treatment, dryrun):
    design = Design(file, control, [control, treatment])
    design.create_design_matrix()
    file_name = directory + '{treatment}_vs_{control}_design_matrix.txt'.format(control=control, treatment=treatment)
    if not dryrun:
        design.design_matrix.to_csv(file_name, index=False, sep="\t")



if __name__ == '__main__':
    readSamples()
