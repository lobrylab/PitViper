import pandas as pd
import sys


def get_conditions(table):
    "Create a dictionary of conditions and column names."
    conditions = {}
    for col in table.columns[2:]:
        if col[:-2] not in conditions.keys():
            conditions[str(col[:-2])] = []
            conditions[str(col[:-2])].append(col)
        elif col[:-2] in conditions.keys():
            conditions[str(col[:-2])].append(col)
    return conditions


def calculate_means(table, baseline):
    "Add a column with mean of control's values."
    conditions = get_conditions(table)
    table[baseline] = table[conditions[baseline]].mean(axis = 1)
    return table


baseline = sys.argv[2]
print("Baseline: ", baseline)


count_table = sys.argv[1]
print("Count table: ", count_table)

output = './data/counts_for_BAGEL.txt'
print("Output: ", output)


data = pd.read_csv(count_table, sep = "\t")

conditions = get_conditions(data)
print(conditions)

data_mean = calculate_means(data, baseline)

cols= conditions[baseline]
data_mean_by_condition = data_mean.drop(cols, axis=1)

print("Writting output...")
data_mean_by_condition.to_csv(output, index = False, sep = "\t")
