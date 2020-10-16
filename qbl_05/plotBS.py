import matplotlib.pyplot as plt

### This is the list of file names, they were counted in Vim
# ER4_exon_overlap
# ER4_intron_overlap
# ER4_promoter_overlap
# G1E_exon_overlap
# G1E_intron_overlap
# G1E_promoter_overlap

er4_overlap = [101,341,67]
g1e_overlap = [85,311,55]
key = ["Exon", "Intron","Promoter"]

sites = [598, 664]
key2 = ["Sites Lost", "Sites Gained"]

fig, ax = plt.subplots(2,1,figsize = (14,10))
ax[0].bar(key,er4_overlap,label='ER4')
ax[0].bar(key,g1e_overlap,label='G1E')
ax[1].bar(key2,sites)
ax[0].set_title('CTCF binding sites')
ax[0].set_xlabel('CTCF')
ax[0].legend()
ax[1].set_title("Number of sites gained/lost")
ax[1].set_ylabel("Number of sites")
plt.savefig("test.png")
