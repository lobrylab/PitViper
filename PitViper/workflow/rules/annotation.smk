rule awk_bed_formating:
	input:
		config['bed_annotation_file']
	output:
		"resources/" + config['token'] + "/annotation_ROSE.bed"
	log:
		"logs/" + config['token'] + "/awk_formating.log"
	shell:
  		"bash ./workflow/scripts/awk_formating.sh {input} {output}"


rule ROSE_annotation:
	input:
		"resources/" + config['token'] + "/annotation_ROSE.bed"
	output:
		"resources/" + config['token'] + "/annotation_ROSE_REGION_TO_GENE.txt"
	params:
		out_dir = "resources/" + config['token'] + "/",
		genome_version = config['genome_version_rose']
	log:
		"logs/" + config['token'] + "/ROSE_annotation.log"
	shell:
		"python workflow/scripts/ROSE/ROSE_geneMapper.py -i {input} -g {params.genome_version} -o {params.out_dir}"