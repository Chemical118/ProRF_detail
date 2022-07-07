from fpbase import *
from selenium import webdriver
from fake_useragent import UserAgent
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import math

root = "avgfp"

ua = UserAgent()
options = webdriver.ChromeOptions()
options.add_argument(f'user-agent={ua.random}')
options.add_argument('headless')
options.add_argument('window-size=1920x1080')
options.add_argument("disable-gpu")
driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)

leaf_list = get_tree_fpbase(driver, root)

ans_list = []
seq_list = []
for name, link in leaf_list[:5]:
    print(name)
    table_list = crawl_dataset(link)
    if len(table_list) > 1:
        raw_list = table_list[1].values.tolist()
        seq, ref = crawl_data(driver, link)
        if len(raw_list) > 1:
            ans_list.append([name, *min(raw_list, key=lambda t: t.count(math.nan))[1:], ref])
        else:
            ans_list.append([name, *raw_list[0], ref])

        seq_list.append(SeqRecord(seq=Seq(seq), id=name, name="", description=""))

driver.quit()

print("Length : %d" % len(ans_list))
print("------")
for li in ans_list:
    print(*li, sep='|')

SeqIO.write(seq_list, "%s_data.fasta" % root, "fasta")
