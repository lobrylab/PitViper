# rule awk_bed_formating:
# 	"""This rule is used to format the annotation file to be used by ROSE."""
# 	input:
# 		config['bed_annotation_file']
# 	output:
# 		f"resources/{config['token']}/annotation_ROSE.bed"
# 	log:
# 		f"logs/{config['token']}/awk_formating.log"
# 	message:
# 		"Formating annotation file {input} to {output}."
# 	shell:
#   		"bash ./workflow/scripts/awk_formating.sh {input} {output}"


rule awk_bed_formating:
	"""This rule is used to format the annotation file to be used by ROSE."""
	input:
		config['bed_annotation_file']
	output:
		f"resources/{config['token']}/annotation_ROSE.bed"
	log:
		f"logs/{config['token']}/awk_formating.log"
	message:
		"Formating annotation file {input} to {output}."
	script:
  		"../../workflow/scripts/bed_to_rose.py"

rule ROSE_annotation:
	"""This rule is used to annotate the regions with ROSE."""
	input:
		f"resources/{config['token']}/annotation_ROSE.bed"
	output:
		f"resources/{config['token']}/annotation_ROSE_REGION_TO_GENE.txt"
	params:
		out_dir = f"resources/{config['token']}/",
		genome_version = config['genome_version_rose']
	log:
		f"logs/{config['token']}/ROSE_annotation.log"
	message:
		"Annotating regions with ROSE: {input} to {output} with genome version {params.genome_version}."
	shell:
		"python workflow/scripts/ROSE/ROSE_geneMapper.py -i {input} -g {params.genome_version} -o {params.out_dir}"
