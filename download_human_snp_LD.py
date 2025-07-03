import sys,os
import requests, time


#snp_list = [x.rstrip() for x in open('hg38_union_loop_leadSNP_100kwindow.txt')]
snp_list = [x.rstrip() for x in open(sys.argv[1])]

token = "720c82118569"  # get from https://ldlink.nci.nih.gov/?tab=apiaccess

for snp in snp_list:
    if os.path.exists('LD_snp/%s.txt'%snp):
        continue
    params = {
        "var": snp,
        "pop": "EUR",
        "r2_d": "r2",
        'window': 100000,
        "genome_build": 'grch38',
        "token": token
    }
    url = "https://ldlink.nih.gov/LDlinkRest/ldproxy"
    r = requests.get(url, params=params)
    if r.status_code == 200:
        #print(f"LD results for {snp}:\n")
        out = open('LD_snp/%s.txt'%snp, 'w')
        out.write(r.text)
        out.close()
    else:
        out = open('no_ld_snp.txt', 'a')
        out.write(snp+'\n')
        out.close()
        #print(f"Error for {snp}: {r.status_code}")
    time.sleep(0.5)  # prevent rate limit
