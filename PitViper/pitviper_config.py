import click
import yaml
import os
import os.path

class ConfigPitViper:
    def __init__(self, token, screen_type, start_from, count_table, samples, files, library, mageck_mle_activate, mageck_mle_norm_method, mageck_rra_activate, bagel_activate, bagel_essentials, bagel_nonessentials, crisphiermix_activate):
        self.token = token
        self.screen_type = screen_type
        self.start_from = start_from
        self.count_table = count_table
        self.samples = samples
        self.files = files
        self.library = library
        self.mageck_mle_activate = mageck_mle_activate
        self.mageck_mle_norm_method = mageck_mle_norm_method
        self.mageck_rra_activate = mageck_rra_activate
        self.bagel_activate = bagel_activate
        self.bagel_essentials = bagel_essentials
        self.bagel_nonessentials = bagel_nonessentials
        self.crisphiermix_activate = crisphiermix_activate

    def show_config(self):
        print('------------------------------')
        print('Configuration file metadata')
        print('------------------------------')
        print('Token:\t', self.token)
        print('Screen type:\t', self.screen_type)
        print('Start from:', self.start_from)
        print('Count table:\t', self.count_table)
        print('Samples file:\t', self.samples)
        print('Files:\t', self.files)
        print('Library:\t', self.library)
        print('MAGeCK MLE activate:\t', self.mageck_mle_activate)
        print('MAGeCK MLE normalization method:\t', self.mageck_mle_norm_method)
        print('MAGeCK RRA activate:\t', self.mageck_rra_activate)
        print('BAGEL activate:\t', self.bagel_activate)
        print('BAGEL essential genes:\t', self.bagel_essentials)
        print('BAGEL nonessential genes:\t', self.bagel_nonessentials)
        print('CRISPhieRmix activate:\t', self.crisphiermix_activate)
        print('------------------------------')


    def write_yaml_conf(self):
        file_out_name = 'config/{token}.yaml'.format(token=self.token)
        print('> Writing yaml config file:', file_out_name)
        dict_file = {'token' : self.token,
                     'screen_type' : self.screen_type,
                     'start_from' : self.start_from,
                     'inputs' : {'count_table' : self.count_table,
                                 'samples' : self.samples,
                                 'files' : self.files,
                                 'library' : self.library},
                     'MAGeCK' : {'MLE' : {'activate' : self.mageck_mle_activate,
                                          'norm-method' : self.mageck_mle_norm_method},
                                 'RRA' : {'activate' : self.mageck_rra_activate}},
                     'BAGEL' : {'activate' : self.bagel_activate,
                                'essentials' : self.bagel_essentials,
                                'nonessentials' : self.bagel_nonessentials},
                     'CRISPhieRmix' : {'activate' : self.crisphiermix_activate}}

        with open(file_out_name, 'w') as file:
            documents = yaml.dump(dict_file, file)

def get_config(form_results):
    parameters_cmd = ""
    parameters = ["token", "screen_type", "count_from", "count_table", "samples", "files", "library", "mageck_mle_activate", "mageck_mle_norm_method", "mageck_rra_activate", "bagel_activate", "bagel_essentials", "bagel_nonessentials", "crisphiermix_activate"]
    for parameter in parameters:
        if parameter in form_results:
            value = form_results[parameter]
        else:
            value = ""
        parameters_cmd += "--{param} {value}".format(param=parameter, value=form_results[parameter])
    return parameters_cmd


@click.command()
@click.option('--token', help='Token.', required=True, type=str)
@click.option('--screen_type', help='Screen type (eg. "shRNA", "CRISPRi")).', required=True, type=str)
@click.option('--start_from', help='Count from (eg. "bam", "fastq", "count")).', required=True, type=str)
@click.option('--count_table', help='Count table.', required=True, type=str)
@click.option('--samples', help='Samples file.', required=True, type=str)
@click.option('--files', help='Files file.', required=True, type=str)
@click.option('--library', help='Library file.', required=True, type=str)
@click.option('--mageck_mle_activate', help='MAGeCK MLE activate', required=True, type=str)
@click.option('--mageck_mle_norm_method', help='MAGeCK MLE normalization method.', required=True, type=str)
@click.option('--mageck_rra_activate', help='MAGeCK RRA activate.', required=True, type=str)
@click.option('--bagel_activate', help='BAGEL activate.', required=True, type=str)
@click.option('--bagel_essentials', help='BAGEL essential genes.', required=True, type=str)
@click.option('--bagel_nonessentials', help='BAGEL nonessential genes.', required=True, type=str)
@click.option('--crisphiermix_activate', help='CRISPhieRmix activate.', required=True, type=str)
def main(token, screen_type, start_from, count_table, samples, files, library, mageck_mle_activate, mageck_mle_norm_method, mageck_rra_activate, bagel_activate, bagel_essentials, bagel_nonessentials, crisphiermix_activate):
    config_pitviper = ConfigPitViper(token, screen_type, start_from, count_table, samples, files, library, mageck_mle_activate, mageck_mle_norm_method, mageck_rra_activate, bagel_activate, bagel_essentials, bagel_nonessentials, crisphiermix_activate)
    config_pitviper.show_config()
    yaml_conf_file = 'PitViper/config/{token}.yaml'.format(token=token)
    if os.path.isfile(yaml_conf_file):
        print("Info: Configuration file {configfile} already exist.".format(configfile=yaml_conf_file))
    else:
        config_pitviper.write_yaml_conf()


if __name__ == '__main__':
    main()
