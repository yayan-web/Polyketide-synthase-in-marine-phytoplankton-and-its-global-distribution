#!/usr/bin/env python
# coding: utf-8

# #### 进化树的修改

# In[952]:


from Bio import SeqIO
import pandas as pd
import numpy as np
from mizani.colors import brewer


# In[309]:


# 定义一个函数来读取FASTA文件并转换为DataFrame
def fasta_to_dataframe(fasta_file):
    # 创建一个空的列表，用于存储序列记录的ID和序列
    sequences = []

    # 遍历FASTA文件中的每个序列记录
    for record in SeqIO.parse(fasta_file, "fasta"):
        # 将记录的ID和序列添加到列表中
        sequences.append({
            'id': record.description,
            'sequence': str(record.seq)
        })

    # 创建DataFrame
    df = pd.DataFrame(sequences)

    return df


# In[310]:


fasta_file = '../浮游植物PKS基因的多样性/不同藻类/进化树/artical/KS_2.fasta'
fasta_df = fasta_to_dataframe(fasta_file)

# 显示DataFrame的前几行
print(fasta_df.head())


# In[311]:


fasta_df.columns


# In[312]:


fasta_df[["ID", "Genus", "Species", "Domain number", "other"]] = fasta_df["id"].str.split(" ", n = 4, expand = True)


# In[313]:


fasta_df.head()


# In[314]:


fasta_df = fasta_df.drop("id", axis = 1)


# In[315]:


tax_df = pd.read_csv("tax_inf.csv", header = 0)
fasta_df = pd.merge(tax_df, fasta_df, on = "Genus")


# In[316]:


fasta_df.head()


# In[317]:


len(fasta_df["sequence"])


# In[318]:


len(set(fasta_df["Species"]))


# In[319]:


cyano = fasta_df[fasta_df["Class"] == "Cyanophyceae"]
cyano["Species"] = cyano["Genus"] + " " + cyano["Species"]


# In[320]:


print(len(cyano["ID"]))
print(set(cyano["Species"]))
print(len(set(cyano["Species"])))


# In[321]:


fasta_df["Species"] = fasta_df["Genus"] + " " + fasta_df["Species"]


# In[322]:


# 设置随机种子
np.random.seed(2024)  # 可以根据需要更改种子值

# 按列名分组，并对每个组应用抽样
def custom_sample(group):
    if len(group) < 3:
        return group
    else:
        return group.sample(n=3)


# 根据某列的元素随机挑选每个元素3个
fasta_df_2 = fasta_df.groupby('Species').apply(custom_sample).reset_index(drop=True)

# 打印新数据框
print(fasta_df_2.head())


# In[323]:


fasta_df_3 = fasta_df_2[["Kingdom", "Phylum", "Class", "Species", "ID", "Domain number", "sequence"]]


# In[324]:


fasta_df_3['id'] = fasta_df_3[["Kingdom", "Phylum", "Class", "Species", "ID", "Domain number"]].apply(lambda row: ' '.join(row), 
                                                                                                      axis=1)


# In[325]:


fasta_file = 'KS—selected.fasta'

# 将数据写入FASTA文件
with open(fasta_file, 'w') as file:
    for index, row in fasta_df_3.iterrows():
        file.write('>' + row['id'] + '\n')  # 写入序列ID
        file.write(row['sequence'] + '\n')  # 写入序列


# ### 进化树绘制

# ##### align：ClustalW
# #### IQ-tree

# In[326]:


# 对进化树的源文件进行修改
def remove_strings_from_file(filename, strings_to_remove):
    with open(filename, 'r') as file:
        lines = file.readlines()

    with open(filename, 'w') as file:
        for line in lines:
            for string in strings_to_remove:
                line = line.replace(string, '')
            file.write(line)


# In[329]:


# 使用示例
filename = 'KS—selected1.nwk' # 替换成你的文件路径
strings_to_remove = ['Bacteria Cyanobacteria Cyanophyceae ', 
                     'Chromista Haptophyta Prymnesiophyceae ', 
                     'Chromista Myzozoa Dinophyceae ',
                     'Chromista Ochrophyta Coscinodiscophyceae ',
                     'Chromista Ochrophyta Pelagophyceae ',
                     'Plantae Chlorophyta Mamiellophyceae ',
                     'Plantae Chlorophyta Trebouxiophyceae ',
                     'Protozoa Amoebozoa Discosea ',
                     'Protozoa Choanozoa Choanoflagellatea ',
                     'Viridiplantae Chlorophyta core chlorophytes ',
                     'Viridiplantae Streptophyta Magnoliopsida '
                    ] # 替换成你想要删除的特定字符串

remove_strings_from_file(filename, strings_to_remove)


# #### 网络图绘制

# In[ ]:


# ugene软件对序列进行对其和裁剪


# In[594]:


from Bio import pairwise2
from Bio import SeqIO
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import math
from sklearn.preprocessing import MinMaxScaler
from mizani.colors import brewer


# In[598]:


# 从文件中读取蛋白序列
filename = 'KS—selected.fasta'

protein_sequences = []
with open(filename, 'r') as file:
    for record in SeqIO.parse(file, 'fasta'):
        protein_sequences.append(record)


# In[599]:


len(protein_sequences)


# In[600]:


# 生成node文件和edge文件
edge_file = 'edges_2.csv'  # 替换成你要保存的edge文件路径
node_file = 'nodes_2.csv'  # 替换成你要保存的node文件路径


# In[601]:


def perform_alignment(sequence1, sequence2):
    alignments = pairwise2.align.globalxx(sequence1, sequence2)
    best_alignment = alignments[0]
    alignment_score = best_alignment.score
    return alignment_score


# In[602]:


num_alignments = math.comb(len(protein_sequences), 2)
edge_data = []
count = 0
completed_percent = 0
for i in range(len(protein_sequences)):
    sequence1 = protein_sequences[i]
    for j in range(i+1, len(protein_sequences)):
        sequence2 = protein_sequences[j]
        alignment_score = perform_alignment(sequence1.seq, sequence2.seq)  # 执行比对，你需要根据具体的比对算法来实现
        edge_data.append((sequence1.description, sequence2.description, alignment_score))
        count += 1
        percent_completed = (count / num_alignments) * 100
        if percent_completed // 5 > completed_percent:
            completed_percent = percent_completed // 5
            print(f"Completed {percent_completed}% of alignments.")


# In[603]:


num_alignments


# In[604]:


# 将edge数据保存到csv文件
df_edges = pd.DataFrame(edge_data, columns=['ID1', 'ID2', 'Alignment Score'])


# In[605]:


node_data = [(sequence.description,) for sequence in protein_sequences]
node_data = pd.DataFrame(node_data)
node_data = node_data.rename(columns={0: 'id'})
node_data[["Kingdom", "Phylum", "Class", "Genus", "Species", "ID", "Domain number"]] = node_data["id"].str.split(" ", n = 6, expand = True)


# In[606]:


node_data.head()


# In[667]:


df_edges.head()


# In[608]:


# 查看Alignment Score的分布
sns.kdeplot(df_edges['Alignment Score'], shade = True)
plt.title('Distribution of Column1')
plt.xlabel('Values')
plt.ylabel('Frequency')
# 显示图形
plt.show()


# In[747]:


# 根据 Alignment Score 筛选edge表格
df_edges_2 = df_edges.copy()


# In[748]:


included_chars = ["Dinophyceae"]
df_edges_2 = df_edges_2[df_edges_2['ID1'].str.contains('|'.join(included_chars)) & df_edges_2['ID2'].str.contains('|'.join(included_chars))]


# In[671]:


# 查看Alignment Score的分布
# sns.kdeplot(df_edges_2['Alignment Score'], shade = True)
plt.hist(df_edges_2['Alignment Score'], bins=50, edgecolor='black')
plt.title('Distribution of Column1')
plt.xlabel('Values')
plt.ylabel('Frequency')
# 显示图形
plt.show()


# In[750]:


df_edges_2 = df_edges_2[df_edges_2["Alignment Score"] > 350]


# In[751]:


node_data_2 = node_data[node_data["Class"] == "Dinophyceae"]


# In[752]:


# 绘制网络图
# 创建空图
G = nx.Graph()

# 添加节点
for i, row in node_data_2.iterrows():
    G.add_node(row['id'], attribute=row['Genus'])

# 添加边
for i, row in df_edges_2.iterrows():
    G.add_edge(row['ID1'], row['ID2'], weight=row['Alignment Score'])


# In[753]:


print("Number of nodes:", G.number_of_nodes())
print("Number of edges:", G.number_of_edges())


# In[754]:


color_palette = {'Alexandrium': '#003666', # 亚历山大藻
                 'Azadinium': '#00aeff',
                 'Gambierdiscus': '#3369e7', # 冈比亚藻
                 'Hematodinium': '#b84592', # 一种寄生甲藻
                 'Karenia': '#ff4f81', # 凯伦藻
                 'Karlodinium': '#ff6c5f', # 卡罗藻属
                 'Prorocentrum': '#ffc168', # 原甲藻属
                 'Symbiodinium': '#2dde98' # 共生甲藻
                }

node_color = [color_palette[G.nodes[n]['attribute']] for n in G.nodes()]


# In[755]:


# 计算中心性
centrality = nx.betweenness_centrality(G, k=10, endpoints=True)
# 计算 community structure
lpc = nx.community.label_propagation_communities(G)
community_index = {n: i for i, com in enumerate(lpc) for n in com}


# In[756]:


# 创建图例
legend_labels = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=10)
                 for label, color in color_palette.items()]
    
plt.legend(handles=legend_labels, title='Node Attribute', loc='best')
plt.savefig('Dinophyceae_ks_net_225_legend.pdf', dpi=300, format='pdf', bbox_inches='tight')


# In[914]:


pos = nx.spring_layout(G, k = 6, seed = 16)
# pos = nx.spectral_layout(G, weight='Alignment Score', scale = 2, dim = 2)
plt.figure(figsize = (4,4))
nx.draw(
    G,
    connectionstyle = "arc3, rad=1",
    pos = pos, 
    edge_color = '#6a737b',
    width = 0.5, 
    with_labels = False,
    node_size = 300, 
    node_color = node_color,
    alpha = 0.8
)


plt.axis('off')
# plt.show()
plt.savefig('Dinophyceae_ks_net_225.tiff', dpi=300, format='tiff', bbox_inches='tight')


# #### 全部点的网络

# 
# ```python
# # 全部的
# df_edges_3 = df_edges.copy()
# 
# 
# # 查看Alignment Score的分布
# # sns.kdeplot(df_edges_2['Alignment Score'], shade = True)
# plt.hist(df_edges_3['Alignment Score'], bins=50, edgecolor='black')
# plt.title('Distribution of Column1')
# plt.xlabel('Values')
# plt.ylabel('Frequency')
# # 显示图形
# plt.show()
# df_edges_3 = df_edges_3[df_edges_3["Alignment Score"] > 350]
# 
# 
# # 绘制网络图
# # 创建空图
# G = nx.Graph()
# 
# # 添加节点
# for i, row in node_data.iterrows():
#     G.add_node(row['id'], attribute=row['Phylum'])
# 
# # 添加边
# for i, row in df_edges_3.iterrows():
#     G.add_edge(row['ID1'], row['ID2'], weight=row['Alignment Score'])
# 
# print("Number of nodes:", G.number_of_nodes())
# print("Number of edges:", G.number_of_edges())
# 
# color_palette = {'Amoebozoa': '#003666', # 变形虫门
#                  'Chlorophyta': '#00aeff', # 绿藻门
#                  'Choanozoa': '#3369e7', # 领鞭毛虫
#                  'Cyanobacteria': '#b84592', # 蓝藻门
#                  'Haptophyta': '#ff4f81', # 囊泡藻界
#                  'Myzozoa': '#ff6c5f', # 粘孢子总门
#                  'Ochrophyta': '#ffc168', # 褐藻门
#                  'Streptophyta': '#2dde98' # 链型植物
#                 }
# 
# node_color = [color_palette[G.nodes[n]['attribute']] for n in G.nodes()]
# 
# 
# # 创建图例
# legend_labels = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=10)
#                  for label, color in color_palette.items()]
#     
# plt.legend(handles=legend_labels, title='Node Attribute', loc='best')
# plt.savefig('Dinophyceae_ks_net_225_legend.pdf', dpi=300, format='pdf', bbox_inches='tight')
# 
# 
# pos = nx.spring_layout(G, k = 7, seed = 2024)
# # pos = nx.spectral_layout(G, weight='Alignment Score', scale = 2, dim = 2)
# plt.figure(figsize = (4,4))
# nx.draw(
#     G,
#     connectionstyle = "arc3,rad=3",
#     pos = pos, 
#     edge_color = '#6a737b',
#     width = 1, 
#     with_labels = False,
#     node_size = 100, 
#     node_color = node_color,
#     alpha = 1
# )
# 
# 
# plt.axis('off')
# plt.show()
# # plt.savefig('all_ks_net_225.tiff', dpi=300, format='tiff', bbox_inches='tight')

# In[917]:


## 删除游离的点


# In[700]:


df_edges_2.head()


# In[ ]:





# #### 根据图总结出两个clade

# In[892]:


components = sorted(nx.connected_components(G), key=len, reverse=True)


# In[ ]:


disconnected_subsets = [component for component in components if len(component) > 1]


# In[896]:


len(disconnected_subsets)


# In[898]:


subset1 = disconnected_subsets[1] # 11
print(len(subset1))
subset2 = disconnected_subsets[2] # 5 Gambierdiscus
print(len(subset2))
subset3 = disconnected_subsets[3] # 2 Karenia
print(len(subset3))
subset4 = disconnected_subsets[0] # 22
print(len(subset4))


# In[924]:


set1 = fasta_df_3[fasta_df_3["id"].isin(subset1)]
print(len(set1))
set2 = fasta_df_3[fasta_df_3["id"].isin(subset2)]
print(len(set2))
set3 = fasta_df_3[fasta_df_3["id"].isin(subset3)]
print(len(set3))
set4 = fasta_df_3[fasta_df_3["id"].isin(subset4)]
print(len(set4))


# In[ ]:


# 将4个文件写出为fasta文件
def write_fasta_file(dataframe, id_column, sequence_column, filename):
    with open(filename, 'w') as file:
        for index, row in dataframe.iterrows():
            sequence_id = row[id_column]
            sequence = row[sequence_column]
            
            file.write(f">{sequence_id}\n")
            file.write(f"{sequence}\n")

write_fasta_file(set1, 'id', 'sequence', 'set1.fasta')
write_fasta_file(set2, 'id', 'sequence', 'set2.fasta')
write_fasta_file(set3, 'id', 'sequence', 'set3.fasta')
write_fasta_file(set4, 'id', 'sequence', 'set4.fasta')


# In[ ]:


# 结构域的数目比较


# In[934]:


set1["component"] = "Set 1"
set2["component"] = "Set 2"
set3["component"] = "Set 3"
set4["component"] = "Set 4"


# In[943]:


set3.head()


# In[938]:


ks_sets = pd.concat([set1, set2, set3, set4], axis = 0)


# In[940]:


ks_sets.head()


# In[997]:


ks_sets["Domain number"] = pd.to_numeric(ks_sets["Domain number"])
p = (ggplot(ks_sets, 
            aes(x = "component", y = "Domain number", color = "component")) +
    scale_color_manual(values=['#cf8d2e', '#2c9f45', '#52325d', "#be0027"]) +
    geom_jitter(shape =  ".", size = 5) +
    # geom_violin() +
    theme_bw() +
    labs(x = "Network component", y = "Domain number") +
    theme(legend_position = "none",
         axis_text = element_text(family = "arial", size = 10),
         axis_title = element_text(family = "arial", size = 15),
         figure_size = (5, 3))
    ) 
p


# In[996]:


ggsave(p, "set_domian_compare.tiff", dpi = 300)


# ### 不同结构域和甲藻的关系

# In[989]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set the working directory
import os
os.chdir("D:\论文\PKS\浮游植物PKS基因的多样性\全球分布")

# Read phy_ann.csv
phy_ann = pd.read_csv("phy_ann.csv", header=0, encoding="UTF-8")
print(phy_ann.columns)
print(phy_ann["Phylum"].unique())

# Subset phy_ann for "Dinophyta" Phylum
dino_ann = phy_ann[phy_ann["Phylum"] == "Dinophyta"][["Phylum", "Order"]].drop_duplicates()

# Read phy_order_data.csv
phy_data = pd.read_csv("phy_order_data.csv", 
                       header=0, encoding="UTF-8")
phy_data = phy_data.rename(columns={phy_data.columns[1]: "Order"})

# Calculate relative abundance
phy_data2 = phy_data.iloc[:, 2:337].transpose()
phy_data2_sum = phy_data2.sum(axis=1)
phy_data3 = phy_data2.divide(phy_data2_sum, axis=0)
phy_data4 = phy_data3.transpose()
phy_data5 = pd.concat([phy_data[["Order"]], phy_data4], axis=1)

# Merge dino_ann and phy_data5
dino_data = pd.merge(dino_ann, phy_data5, on="Order")
dino_data = dino_data.transpose()
dino_data.columns = dino_data.iloc[0]
dino_data = dino_data[1:]


# In[990]:


# Read domin_env_reads_2.csv
prd_data = pd.read_csv("domin_env_reads_2.csv", header=0, encoding="UTF-8")


# In[994]:


prd_data.columns[11] = "Barcode"


# In[988]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set the working directory
import os
os.chdir("D:\论文\PKS\浮游植物PKS基因的多样性\全球分布")

# Read phy_ann.csv
phy_ann = pd.read_csv("phy_ann.csv", header=0, encoding="UTF-8")
print(phy_ann.columns)
print(phy_ann["Phylum"].unique())

# Subset phy_ann for "Dinophyta" Phylum
dino_ann = phy_ann[phy_ann["Phylum"] == "Dinophyta"][["Phylum", "Order"]].drop_duplicates()

# Read phy_order_data.csv
phy_data = pd.read_csv("phy_order_data.csv", 
                       header=0, encoding="UTF-8")
phy_data = phy_data.rename(columns={phy_data.columns[1]: "Order"})

# Calculate relative abundance
phy_data2 = phy_data.iloc[:, 2:337].transpose()
phy_data2_sum = phy_data2.sum(axis=1)
phy_data3 = phy_data2.divide(phy_data2_sum, axis=0)
phy_data4 = phy_data3.transpose()
phy_data5 = pd.concat([phy_data[["Order"]], phy_data4], axis=1)

# Merge dino_ann and phy_data5
dino_data = pd.merge(dino_ann, phy_data5, on="Order")
dino_data = dino_data.transpose()
dino_data.columns = dino_data.iloc[0]
dino_data = dino_data[1:]

# Read domin_env_reads_2.csv
prd_data = pd.read_csv("domin_env_reads_2.csv", header=0, encoding="UTF-8")
prd_data = prd_data.rename(columns={prd_data.columns[11]: "Barcode"})

# Merge sample_ids and dino_data
sample_ids = pd.read_csv("enviro_18SV9v1.csv", sep=",", header=0)
dino_data = pd.merge(sample_ids[["Barcode", "Stations"]], dino_data, on="Barcode")

# Merge prd_data and dino_data
prd_env_dino = pd.merge(prd_data, dino_data, on="Stations")

# Write prd_env_dino to CSV
prd_env_dino.to_csv("prd_env_dino.csv", index=False)

# Plot using ggplot
import plotnine as pn
from plotnine import *
from plotnine.data import *

plt.figure(figsize=(10, 6))
p = (
    ggplot(prd_env_dino_long)
    + geom_point(aes(x="Gonyaulacales", y="value", color="variable"))
    + facet_wrap("~variable", scales="free_y")
)
print(p)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ### 在linux服务器上运行

# In[ ]:





# In[844]:


# ugene对序列进行对齐
# hmmbuild output.hmm input.fasta # in linux shell

# idx name                  nseq  alen  mlen eff_nseq re/pos description
#---- -------------------- ----- ----- ----- -------- ------ -----------
1     set1                    11   240   226     0.99  0.590 

# idx name                  nseq  alen  mlen eff_nseq re/pos description
#---- -------------------- ----- ----- ----- -------- ------ -----------
1     set2                     5   383   372     0.59  0.589 

# idx name                  nseq  alen  mlen eff_nseq re/pos description
#---- -------------------- ----- ----- ----- -------- ------ -----------
1     set3                     2   390   386     0.55  0.589 

# idx name                  nseq  alen  mlen eff_nseq re/pos description
#---- -------------------- ----- ----- ----- -------- ------ -----------
1     set4                    22   429   376     1.74  0.590 


# In[ ]:


# hmmsearch --F3 1e-3 -o set1.txt set1.hmm ../tara_data/trans/seqs/MATOU-v1.fna
# hmmsearch --F3 1e-3  -o set4.txt set4.hmm ../tara_data/trans/seqs/MATOU-v1.fna


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[289]:


import pandas as pd
from plotnine import *
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import warnings

warnings.filterwarnings("ignore")


# In[290]:


net_df = pd.read_csv("KSs Full Network_1 default node.csv", header = 0)
tax_df = pd.read_csv("tax_inf.csv", header = 0)


# In[291]:


net_df = net_df[["shared name", "Description", "Sequence Status", "Sequence Length"]]
new_df = net_df["Description"].str.split(expand = True)
new_df.rename(columns = {0:"ID", 1:"Genus", 2:"Species", 3:"Domain number", 4:"Other"},
             inplace = True)


# In[292]:


df = pd.concat([net_df, new_df], axis = 1)
set(df["Genus"]) - set(tax_df["Genus"]) # set()
df = pd.merge(df, tax_df, on = "Genus", how = "left")
df.to_csv('node_inf.csv', index=False)
# df.columns
new_column_order = ['shared name', 'Desscription', 'Sequence Status', "Sequence Length",
                   "Domain number", "Kingdom", "Phylum", "Class", "Genus", "Species"]

df_2= df.reindex(columns=new_column_order)


# In[293]:


df_2['Domain number'] = df_2['Domain number'].replace('domains', pd.NA)
df_2["Domain number"] = pd.to_numeric(df_2["Domain number"])


# In[294]:


df_2.columns


# In[295]:


df_3 = df_2[df_2["Sequence Status"] == "complete"]


# In[296]:


pairwise_tukey = pairwise_tukeyhsd(df_3['Sequence Length'], df_3['Kingdom'])
print(pairwise_tukey)


# In[297]:


p = (ggplot(df_3, aes(x = "Kingdom", y = "Sequence Length", fill = "Kingdom")) +
     geom_jitter(shape = ".") +
     geom_boxplot(alpha = 0.6) +
     theme(figure_size = (6, 4)) +
     theme_bw() +
     ylim(0, 500) +
     theme(axis_text_x = element_text(angle = 45, vjust = 1, hjust = 1)))
p


# ### 绘制网络图

# In[298]:


import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns


# In[299]:


node_df = pd.read_csv("node_df_complete.csv", header = 0)
edge_df = pd.read_csv("edge_inf_complete.csv", header = 0)


# In[300]:


node_df.head(2)


# In[301]:


edge_df[["source", "target"]] = edge_df["name"].str.split(",", expand = True)


# In[302]:


edge_df.head(2)


# In[303]:


# 创建空图
G = nx.Graph()

# 添加节点
for i, row in node_df.iterrows():
    G.add_node(row['shared name'], attribute=row['Phylum'])

# 添加边
for i, row in edge_df.iterrows():
    G.add_edge(row['source'], row['target'], weight=row['alignment_score'])


# In[304]:


print("Number of nodes:", G.number_of_nodes())
print("Number of edges:", G.number_of_edges())


# In[305]:


# 简化网络

# num_to_remove = int(len(G) / 1.5)
# nodes = sample(list(G.nodes), num_to_remove)
# G.remove_nodes_from(nodes)

# remove low-degree nodes
# low_degree = [n for n, d in G.degree() if d < 10]
# G.remove_nodes_from(low_degree)

# largest connected component
# components = nx.connected_components(G)
# largest_component = max(components, key=len)
# H = G.subgraph(largest_component)


# In[306]:


# 删除具有特定属性的节点
nodes_to_remove = [node for node, 
                   data in G.nodes(data=True) if 'attribute' in data and data['attribute'] in ["Amoebozoa",
                                                                                         "Choanozoa",
                                                                                         "Myzozoa",
                                                                                         "Streptophyta"]]
G.remove_nodes_from(nodes_to_remove)

# 删除具有特定属性的边
# edges_to_remove = [(u, v) for u, v, data in G.edges(data=True) if 'weight' in edge_df and edge_df['weight'] > 0.5]
# G.remove_edges_from(edges_to_remove)


# In[307]:


# 计算中心性
centrality = nx.betweenness_centrality(G, k=10, endpoints=True)
# 计算 community structure
lpc = nx.community.label_propagation_communities(H)
community_index = {n: i for i, com in enumerate(lpc) for n in com}


# In[ ]:


atts = [G.nodes[n]['attribute'] for n in G.nodes()]
set(atts)


# In[ ]:


color_palette = {# 'Amoebozoa': '#003666', # 变形虫门
                 'Chlorophyta': '#00aeff', # 绿藻门
                 # 'Choanozoa': '#3369e7', # 领鞭毛虫
                 'Cyanobacteria': '#b84592', # 蓝藻门
                 'Haptophyta': '#ff4f81', # 囊泡藻界
                 # 'Myzozoa': '#ff6c5f', # 粘孢子总门
                 'Ochrophyta': '#ffc168', # 褐藻门
                 # 'Streptophyta': '#2dde98' # 链型植物
                }

node_color = [color_palette[G.nodes[n]['attribute']] for n in G.nodes()]


# In[ ]:


# edge_widths = [edge_df['alignment_score'] for u, v, edata in G.edges(data=True)]


# In[ ]:


# 创建图例
legend_labels = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=10)
                 for label, color in color_palette.items()]
plt.legend(handles=legend_labels, title='Node Attribute', loc='best')

pos = nx.spring_layout(G, k = 3, seed = 1)
plt.figure()
nx.draw(
    G, 
    pos = pos, 
    edge_color = '#6a737b',
    width = 0.05, 
    with_labels = False,
    node_size = 50, 
    node_color = node_color,
    alpha = 0.4
)

plt.axis('off')
plt.show()


# In[ ]:




