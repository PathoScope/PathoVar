import re
import itertools
from collections import defaultdict

def scores_as_table(entry_list, top_n = 20):
    entry_list.sort(key = lambda x: x['composite_score'], reverse = True)
    headers = [i for i in itertools.chain.from_iterable([['title'], ['composite_score'], ['snp_burden'], 
        entry_list[0]['blast_scores'].keys(), ["coverage_score"]])]
    lines =  [ [ [s['title']], [s['composite_score']], [s['snp_burden']], 
        s['blast_scores'].values(), [s['coverage_score']]] for s in entry_list[:top_n] ]
    top_entries = [headers]
    for line in lines:
        entry = [i for i in itertools.chain.from_iterable(line)]
        top_entries.append(entry)
    return top_entries

def write_table(entry_list, file_name):
    import csv
    scores_table = scores_as_table(entry_list)
    with open(file_name, 'w') as f:
        w = csv.writer(f, delimiter ='\t')
        w.writerows(scores_table)

def test_heuristic(annotation_report, file_name, blast_max = 20, snp_max = 20, coverage_max = 20,
                        var_score_dict=None, blast_value = 2, coverage_value = 1):
    annotation_report.score_all_entries(blast_max, snp_max, coverage_max,
                        var_score_dict, blast_value, coverage_value)
    sdata = sorted([g for g in annotation_report], key=lambda x: x['composite_score'], reverse=True)
    write_table(sdata, "test_heuristic-snp-%d-blast-%d-%d-cov-%d-%d.tsv" % (snp_max, blast_value, blast_max, coverage_value, coverage_max))




# if missense moderate, 1, if high, 2, based on snpEff? let frame_shift be just high (2)?
# Assume UNKNOWN is missense
# Add a limit on how many mutations can contribute?
# Split score into variant_score ("snp_burden") and blast hits ("drug burden")
# Normalize snp_burden by sequence size
# Add tapering for redundant blast hits.
# split drugbank and card # done
# Max all by 20 and don't normalize by gene size?
# Composite score gets extra multiplier?
# Add weight by allele frequency?
def score_heuristic_cap(entry, blast_max = 20, snp_max = 20, coverage_max = 20,
                        var_score_dict=None, blast_value = 2, coverage_value = 1):
    if var_score_dict is None:
        var_score_dict = {"LOW": 0.1, "UNKNOWN": 1,"MODERATE": 1, "MODIFIER": 1, "HIGH": 2}
    variant_score = 0
    for variant in entry["variants"]:
        variant_score += var_score_dict[variant["eff"].get("impact", "UNKNOWN")]

    entry['snp_burden'] = variant_score
    composite = max(min(variant_score, snp_max), 1)

    blast_db_scores = defaultdict(int)
    current_db = "UNKNOWN"
    for blast_db, blast_hits in entry['blast_hits'].items():
        if re.search("drugbank", blast_db):
            current_db = 'drugbank'
        elif re.search("comprehensive_antibiotic", blast_db):
            current_db = "card"
        else:
            current_db = "UNKNOWN"
        for hit in blast_hits:
            blast_db_scores[current_db + "_score"] += 1

    entry['blast_scores'] = blast_db_scores
    
    blast_total = 0
    for blast_db, score in blast_db_scores.items():
        blast_total += min(score * blast_value, blast_max)

    composite *= max(blast_total, 1)

    coverage_score = 0
    for region in entry['uncovered_regions']:
        coverage_score += coverage_value

    entry['coverage_score'] = coverage_score

    composite *= max(min(coverage_score, coverage_max), 1)

    entry['composite_score'] = composite
    return entry



def overfit_score_heuristic(entry, snp_score = .1, missense_score = 2, frame_shift_score = 2/3.0, min_quality = 20, \
    blast_hit_score = 2, multi_drug_bonus = 2, uncov_region_score = 2):
    variant_score = 1
    seq_len = float(entry['end'] - entry['start'])
    for variant in entry["variants"]:
        if variant['call_quality'] < min_quality:
            continue
        if variant["eff"]["type"] in ["UNKNOWN", "SYNONYMOUS_CODING", "INTERGENIC"]:
            variant_score += snp_score

        elif variant["eff"]["type"] in ["NON_SYNONYMOUS_CODING", "NON_SYNONYMOUS_START", \
                                        "SYNONYMOUS_STOP", "NON_SYNONYMOUS_STOP"]:
            variant_score += missense_score

        elif variant["eff"]["type"] in ["STOP_GAINED", "STOP_LOST", "START_LOST", "RARE_AMINO_ACID"]:
            variant_score += seq_len/3.0 * min(seq_len / 2000, 1)

        elif variant["eff"]["type"] in ["FRAME_SHIFT", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR"]:
            start_point = variant['start'] - entry['start']
            percent_shift = (start_point / seq_len)
            # The earlier the shift occurs, the greater the score multiplier, 
            # but penalize short sequences ( < 1000 bp) 
            score = (frame_shift_score / percent_shift) * min(seq_len / 2000, 1)
            variant_score += score
            #variant_score += frame_shift_score
        else:
            print("Variant Effect Not Recognized: %s" % variant["eff"]["type"])
    
    blast_score = 1
    for blast_db, blast_hits in entry['blast_hits'].items():
        for hit in blast_hits:
            blast_score += blast_hit_score
            if re.search(r'multidrug', hit['hit_def'], re.IGNORECASE):
                blast_score += multi_drug_bonus

    coverage_score = 1
    for region in entry["uncovered_regions"]:
        # if(entry['mean_coverage'] == 0):
        #     uncov_region_score * --
        coverage_score += uncov_region_score

    entry["score"] = variant_score * blast_score * coverage_score
    return entry