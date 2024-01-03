import pandas as pd 
import os 
os.chdir("/hd2/patrick/PXD020394_SARS_COVID19")
os.listdir()
df = pd.read_csv("sdrf_cleaned.csv", sep = "\t")
data_dir = "/hd2/patrick/PXD020394_SARS_COVID19/ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD020394/"
df["comment[file uri]"] = data_dir + df["comment[file uri]"].map(lambda x:x.split("/")[-1]) 


#df = df.rename({"comment[file uri]":"comment[data file]"}, axis = 1)                           
df.to_csv("sdrf_cleaned_local.csv", sep = "\t", index = False)

