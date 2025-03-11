# python3 vcf_compare_location.py -b base.vcf -bn colname_in_b -c call.vcf -cn colname_in_c -gt false -d 0
import argparse
import sys
import os


def read_vcf(loc, coln, consider_genotype, ignore_chr=False):
	loc_dict = dict()
	used_col = 0
	with open(loc) as f:
		for line in f:
			if line.startswith('#'):
				if line.startswith('#CHROM'):
					#temp = line.rstrip().split('\t')
					#for col in range(9, len(temp)):
					#	if colname_include in temp[col]:
					used_col = coln - 1
							#print('#' + loc + ' ' + str(colname_include) + ' ' + str(used_col))
				continue
			else:
				temp = line.rstrip().split('\t')
				if consider_genotype == 1:
					if temp[used_col] in ('.', '0'):
						continue
					if temp[used_col].startswith('./.'):
						continue
					if temp[used_col].startswith('0/0'):
						continue
					if temp[used_col].startswith('0|0'):
						continue
				if ignore_chr:
					chrn = '_'
				else:
					chrn = temp[0]
				pos = int(temp[1])
				if chrn not in loc_dict:
					loc_dict[chrn] = dict()
				if pos not in loc_dict[chrn]:
					loc_dict[chrn][pos] = temp[used_col].split(":")[0]
	return loc_dict

def count_n(dict_set):
	count = 0
	for key in dict_set:
		count += len(list(dict_set[key].keys()))
	return count

def count_redundant_n(dict_set, max_dis):
	count = 0
	for key in dict_set:
		last_pos = -1000000
		for j in dict_set[key]:
			if j - last_pos <= max_dis:
				count += 1
			last_pos = j
	return count

def compare_loc_dict(locb, locc, max_dist):
	tp = 0
	tp_pos = dict()
	for keyb in locb:
		if keyb not in tp_pos:
			tp_pos[keyb] = dict()
		if keyb in locc:
			if max_dist == 0:
				intersect_pos = list(set(locb[keyb].keys()) & set(locc[keyb].keys()))
				intersect_n = len(intersect_pos)
				tp += intersect_n
				for pos in intersect_pos:
					tp_pos[keyb][pos] = ''
			else:
				for posb in locb[keyb]:
					for posc in locc[keyb]:
						if abs(posb - posc) <= max_dist:
							if posb not in tp_pos[keyb] and posc not in tp_pos[keyb]:
								tp += 1
							tp_pos[keyb][posb] = ''
							tp_pos[keyb][posc] = ''
						if posc - posb > max_dist:
							break
	locb_n = count_n(locb) - count_redundant_n(locb, max_dist)
	locc_n = count_n(locc) - count_redundant_n(locc, max_dist)
	fp = locc_n - tp
	tn = locb_n - tp
	precision = 0
	recall = 0
	f1 = 0
	if tp + fp != 0:
		precision = tp / (tp + fp)
	if tp + tn != 0:
		recall = tp / (tp + tn)
	if recall + precision != 0:
		f1 = 2 * ((recall * precision) / (recall + precision))
	jaccard = tp / (tp + fp + tn)
	#print(locb, locc)
	return '\t'.join([str(x) for x in [locb_n, locc_n, tp, fp, tn, round(precision,4), round(recall,4), round(f1,4), round(jaccard,4)]]), tp_pos

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='''Compare vcf location with or without genotype''')
	parser.add_argument('-b', '--base_vcf', metavar='<base.vcf>', help='Path of base.vcf', type=str,
                        required=True)
	parser.add_argument('-bn', '--base_vcf_coln', metavar='<int>', help='Colnumber of base.vcf (default: 10)',
                        type=int, default=10)
	parser.add_argument('-c', '--call_vcf', metavar='<base.vcf>', help='Path of call.vcf', type=str,
						required=True)
	parser.add_argument('-cn', '--call_vcf_coln', metavar='<int>', help='Colnumber of call.vcf (default: 10)',
						type=int, default=10)
	parser.add_argument('-gt', '--consider_genotype', help='Consider genotype (default: 1=Yes)',
						type=int, default=1)
	parser.add_argument('-md', '--max_dis', metavar='<int>', help='Max distance (default: 0)', type=int, default=0)
	parser.add_argument('-o', '--output_vcf', metavar='<output.diff.vcf>', help='Path of output diff vcf [Call-Base] (default: None)', type=str, default='')
    #parser.add_argument('-n', '--sample_name', metavar='<str>', help='Name of sample', type=str, required=True)
	args = vars(parser.parse_args())

	ignore_chr = False

	print('# Input: ' + str(args))
	vcfb = read_vcf(args["base_vcf"], args["base_vcf_coln"], args["consider_genotype"], ignore_chr)
	# print('# Finish reading base vcf.')
	vcfc = read_vcf(args["call_vcf"], args["call_vcf_coln"], args["consider_genotype"], ignore_chr)
	# print('# Finish reading call vcf.')
	res_str, tp_pos_dict = compare_loc_dict(vcfb, vcfc, max_dist=args["max_dis"])
	print('\t'.join(['#Base_vcf', 'Call_vcf' ,'locb_n', 'locc_n', 'tp', 'fp', 'tn', 'precision', 'recall', 'f1', 'jaccard']))
	print(args["base_vcf"] + '\t' + args["call_vcf"] + '\t' + res_str)

	if args["output_vcf"]:
		if os.path.exists(args["output_vcf"]):
			print('# Warning! Output vcf exists!')

		fout = open(args["output_vcf"],'w')
		with open(args["call_vcf"]) as fin:
			for line in fin:
				if line.startswith('#'):
					fout.write(line)
				else:
					temp = line.rstrip().split('\t')
					chrn = temp[0]
					pos = int(temp[1])
					if chrn in tp_pos_dict:
						if pos in tp_pos_dict[chrn]:
							continue
					fout.write(line)
		fout.close()




