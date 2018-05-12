import sys
import urllib2
import time
import numpy as np
import csv

def get_list_of_interacting_proteins(protein, min_score = 900):
    # return a dictionary, with keys as genes and values as number of interactions
    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "interaction_partners"

    # my_genes = ["9606.ENSP00000000233", "9606.ENSP00000000412",
    #            "9606.ENSP00000000442", "9606.ENSP00000001008"]
    my_genes = [protein]
    species = "9606"
    limit = 5
    required_score = min_score
    my_app = "www.awesome_app.org"

    ## Construct the request

    request_url = string_api_url + "/" + output_format + "/" + method + "?"
    request_url += "identifiers=%s" % "%0d".join(my_genes)
    request_url += "&" + "species=" + species
    # request_url += "&" + "limit=" + str(limit)
    request_url += "&" + "required_score=" + str(required_score)
    request_url += "&" + "caller_identity=" + my_app

    try:
        response = urllib2.urlopen(request_url)
    except urllib2.HTTPError as err:
        error_message = err.read()
        print error_message
        sys.exit()

    ## Read and parse the results
    interacting_proteins = []

    line = response.readline()

    while line:
        l = line.strip().split("\t")
        query_ensp = l[0]
        query_name = l[2]
        partner_ensp = l[1]
        partner_name = l[3]
        combined_score = l[5]

        interacting_proteins.append(partner_name)
        #print line
        #print "\t".join([query_ensp, query_name, partner_name, combined_score])

        line = response.readline()

    return interacting_proteins

def get_number_of_interacting_proteins(proteins, min_score = 900):
    # return a dictionary, with keys as genes and values as number of interactions
    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "interaction_partners"

    # my_genes = ["9606.ENSP00000000233", "9606.ENSP00000000412",
    #            "9606.ENSP00000000442", "9606.ENSP00000001008"]
    my_genes = proteins
    species = "9606"
    limit = 5
    required_score = min_score
    my_app = "www.awesome_app.org"

    ## Construct the request

    request_url = string_api_url + "/" + output_format + "/" + method + "?"
    request_url += "identifiers=%s" % "%0d".join(my_genes)
    request_url += "&" + "species=" + species
    # request_url += "&" + "limit=" + str(limit)
    request_url += "&" + "required_score=" + str(required_score)
    request_url += "&" + "caller_identity=" + my_app

    try:
        response = urllib2.urlopen(request_url)
    except urllib2.HTTPError as err:
        error_message = err.read()
        print error_message
        sys.exit()

    ## Read and parse the results
    ret_dict = dict.fromkeys(my_genes, 0)

    line = response.readline()

    while line:
        l = line.strip().split("\t")
        query_ensp = l[0]
        query_name = l[2]
        partner_ensp = l[1]
        partner_name = l[3]
        combined_score = l[5]

        ret_dict[query_name] += 1
        #print "\t".join([query_ensp, query_name, partner_name, combined_score])

        line = response.readline()

    return ret_dict

with open("data/Homo Sapiens - all genes.txt") as f:
    all_human_genes = f.readlines()
all_human_genes = [x.strip() for x in all_human_genes]
all_human_genes = set(all_human_genes)

min_score = 900

with open('human_genome_interactions.csv', 'wb') as csvfile:
    fieldnames = ['Gene', 'Degree (min_score = %d)' % min_score, 'Avg. deg. of neighbors', 'Median deg. of neighbors']
    writer = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(fieldnames)
    for gene in all_human_genes:
        print gene
        try:
            interacting_proteins = get_list_of_interacting_proteins(gene, min_score)
            if not interacting_proteins:
                writer.writerow([gene, 0, 0, 0])
                csvfile.flush()
            else:
                neighbor_of_interacting_proteins = get_number_of_interacting_proteins(
                    interacting_proteins, min_score)
                writer.writerow([gene,
                                len(interacting_proteins),
                                "%.2f" % np.mean(neighbor_of_interacting_proteins.values()),
                                np.median(neighbor_of_interacting_proteins.values())]                            )
                csvfile.flush()
        except:
            print "Error on %s" % gene
            continue
        time.sleep(1) # being nice to STRING