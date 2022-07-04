import ChromosomeMetrics
import matplotlib.pyplot as plt
import Plotting_Methods
import numpy as np


def moving_average(x, w, i):
	if i>w:
		start=i-w
	else:
		start=0
	if i+w > len(x)-1:
		end=len(x)-1+1
	else:
		end=i+w+1
	
	return np.sum(x[start:end])/len(x[start:end])
	#np.convolve(x, np.ones(w), "valid") / w



def PlotBED (bed_file, axis, color, color2):
	x_start=[]
	x_end=[]
	y=[]
	for line in open(bed_file):
		#print(line)
		#print(line.split("\t"))
		
		splits = line.split("\t")
		
		#print(splits[0])
		
		try:
			Chromosom =splits[0]
			Start =int(splits[1])
			Ende=int(splits[2])
			Basen=float(splits[3])

			length = Ende - Start
			
		#	print(length)
			Coverage=Basen/length
			y.append(Coverage)
		#	print(Coverage)
			
			absStart=ChromosomeMetrics.GetAbsolutePosition(Chromosom, Start)
			absStop=ChromosomeMetrics.GetAbsolutePosition(Chromosom, Ende)
			
			x_start.append(absStart)
			x_end.append(absStop)
			
			axis.plot([absStart, absStop], [Coverage, Coverage], c=color)
		#	print(absStart)
		#	print(absStop)
		except ValueError:
		#	print(line)
			break

		#print(moving_average(y, 5))
	
	for i in range(len(y)):		
		axis.plot([x_start[i], x_end[i]], [moving_average(y, 5, i), moving_average(y, 5, i)], c=color2)
	


ax1 = plt.subplot(611)
ax2 = plt.subplot(612)
ax3 = plt.subplot(613)
ax4 = plt.subplot(614)
ax5 = plt.subplot(615)
ax6 = plt.subplot(616)

Plotting_Methods.HighlightEvenChromosomes(ax1, 0, 25)
Plotting_Methods.HighlightEvenChromosomes(ax2, 0, 25)
Plotting_Methods.HighlightEvenChromosomes(ax3, 0, 25)
Plotting_Methods.HighlightEvenChromosomes(ax4, 0, 25)
Plotting_Methods.HighlightEvenChromosomes(ax5, 0, 25)
Plotting_Methods.HighlightEvenChromosomes(ax6, 0, 25)

#PlotBED("../Mappings/Run789_790_791_793-C053_CNV.bed", ax, "b")
PlotBED("../Gliom/Gliom-Amplikons/Run822-1-hg38-chromosomband.bed", ax1, "r", "g")
PlotBED("../Gliom/Gliom-Amplikons/Run822-2-hg38-chromosomband.bed", ax2, "black", "g")
PlotBED("../Gliom/Gliom-Amplikons/Run822-3-hg38-chromosomband.bed", ax3, "b", "g")
PlotBED("../Gliom/Gliom-Amplikons/Run822-4-hg38-chromosomband.bed", ax4, "c", "g")
PlotBED("../Gliom/Gliom-Amplikons/Run822-6-hg38-chromosomband.bed", ax5, "firebrick", "g")
PlotBED("../Gliom/Gliom-Amplikons/Run822-7-hg38-chromosomband.bed", ax6, "y", "g")
#PlotBED("../Gliom/wgs/Run990-2-10kB.bed", ax, "r", "g")
#PlotBED("../Gliom/wgs/Run990-3-1kb.bed", ax, "r", "g")
#PlotBED("../WGS_Test_Mapping/Run713_714_789_790_791_793_876-ngmlr-hg38-chromosomband.bed", ax, "g")
#PlotBED("../WGS_Test_Mapping/WGS.longread-hg38-chromosomband.bed", ax, "b")
#PlotBED("../WGS_Test_Mapping/WGS.longread-wgs_1MB_intervals.bed", ax, "m")
#PlotBED("../WGS_Test_Mapping/Run713_714_789_790_791_793_876-C053_CNV.bed", ax, "b")
#PlotBED("../WGS_Test_Mapping/Run876-wgs_1MB_intervals.bed", ax, "k")
#PlotBED("../Mappings/hg38-chromosomband-run393(2).bed", ax, "g")

ax1.plot([-1, -2], [1, 2], c="r", label="G-BS")
ax2.plot([-1, -2], [1, 2], c="lime", label="G-LE")
ax3.plot([-1, -2], [1, 2], c="black", label="G-FD")

#plt.plot([-1, -2], [7, 8], c="b", label="C053_CNV")

for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
	ax.set_ylabel("Coverage")
	
	
	ax.set_xlim(0, 3088269832)
ax3.set_xlabel("Genom")
ax1.set_ylim(0, 0.01)
ax2.set_ylim(0, 0.006)
ax3.set_ylim(0, 0.02)
ax4.set_ylim(0, 0.002)
ax5.set_ylim(0, 0.1)
ax6.set_ylim(0, 0.08)
#plt.legend()
#plt.title("Coverage C053")
ax6.set_xlabel("Genom")
plt.show()
