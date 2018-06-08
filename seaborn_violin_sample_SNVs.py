import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd    

plt.rcParams["font.family"] = "serif"

#read in the data file from an excel spreadsheet
for SAG in ['C09','E23','K20','M21','N22']:
	data=pd.read_csv("./SeqObj/"+SAG+"_violin_SAAVs.txt",sep='\t')
	#I will do one for Piccard and one for Von Damm


	#use seaborn to make the violin lot
	sns.set_style("ticks") #no grids, just ticks
	sns.despine(trim=True) #removes the top and side lines in the graph so it's just the x and y axis lines drawn
	sns.set_context("paper", font_scale=1.50) #made it suitable for publishing a paper and also changes the font size
	sns.set_color_codes()
	fig = plt.figure(figsize=(10,10)) #i think this adjusts the size
	ax = fig.add_subplot(1,1,1) #no idea what this does but i think it's related to the line above

	#now actually load the data
	ax = sns.violinplot(y="SNVs", x="orf sag/sample", data=data, scale="width")
	ax = sns.violinplot(y="SNVs", x="orf sag/sample", data=data, scale="width")

	ax.tick_params(labelsize=15) #makes the font smaller on the x axis
	locs, labels = plt.xticks() #this and the next line make the labels rotate 45 degrees
	plt.setp(labels, rotation=30)
	plt.savefig(SAG+'_violin.pdf')


#http://stackoverflow.com/questions/31063005/how-to-widen-boxes-in-seaborn-boxplot
#reason why "hue" didn't work. still confused, and i wound up just doing it in illustrator.
