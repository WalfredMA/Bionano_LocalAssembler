#!/usr/bin/python

import pandas as pd 
import StringIO
import re
import os
import argparse
import matplotlib.pyplot as plt
import collections as cl

global stdshift_penalty,cluster_max,base_penalty

def makereverse(a):
	
	allcoordi=sum(a,())[::-1]
	return [(allcoordi[2*k],allcoordi[2*k+1]) for k in xrange(len(allcoordi)/2)]

def calculatedistance(a1,a2):
	
	if (a1==-1 and a2==-1):
		return (a1-a2)**2
	
	
	if min(a1,a2)==0 or max(a1,a2)==1000000:
		
		return 0
	
	else:
		
		return (a1-a2)**2

def similarity(x,y,strand_return=0):
		
	minisize=2
	
	if min(len(x),len(y))<minisize:
		return base_penalty
	
	
	front_x=list(x[0])
	back_x=list(x[-1])
	
	front_y=list(y[0])
	back_y=list(y[-1])
	
	if x[0][1]>=x[0][0]:
		front_x[0]=0
	else:
		front_x[0]=1000000000
	
	if y[0][1]>=y[0][0]:
		front_y[0]=0
	else:
		front_y[0]=1000000000
	
	if x[-1][1]>=x[-1][0]:
		back_x[1]=1000000000
	else:
		back_x[1]=0
	
	if y[-1][1]>=y[-1][0]:
		back_y[1]=1000000000
	else:
		back_y[1]=0
	
	breaks_x=[tuple(front_x)]+x[1:-1]+[tuple(back_x)]
	breaks_y=[tuple(front_y)]+y[1:-1]+[tuple(back_y)]
	
	turn=0
	indent=0
	reverse_x=makereverse(breaks_x)
	
	if strand_return==0:
		
		#differences=min([sum([ calculatedistance(a[0],b[0]) + calculatedistance(a[1],b[1]) if (a[1]-a[0])*(b[1]-b[0])>0 else stdshift_penalty for a,b in zip(breaks_x[max(-turn,0):],breaks_y[max(turn,0):])]) for turn in  xrange(minisize-len(breaks_x),len(breaks_y)-minisize+1)])
		
		#differences_reverse=min([sum([ calculatedistance(a[0],b[0]) + calculatedistance(a[1],b[1]) if (a[1]-a[0])*(b[1]-b[0])>0 else stdshift_penalty for a,b in zip(reverse_x[max(-turn,0):],breaks_y[max(turn,0):])]) for turn in  xrange(minisize-len(reverse_x),len(breaks_y)-minisize+1)])
		
		differences=sum([ calculatedistance(a[0],b[0]) + calculatedistance(a[1],b[1]) if (a[1]-a[0])*(b[1]-b[0])>0 else calculatedistance(a[1],b[0]) + calculatedistance(a[0],b[1])+stdshift_penalty for a,b in zip(breaks_x,breaks_y)]) 
		
		differences_reverse=sum([ calculatedistance(a[0],b[0]) + calculatedistance(a[1],b[1]) if (a[1]-a[0])*(b[1]-b[0])>0 else calculatedistance(a[1],b[0]) + calculatedistance(a[0],b[1])+stdshift_penalty for a,b in zip(reverse_x,breaks_y)]) 
		
				
		return min(differences,differences_reverse)
		
	else:
		
		"""
		#differences=[sum([ calculatedistance(a[0],b[0]) + calculatedistance(a[1],b[1]) if (a[1]-a[0])*(b[1]-b[0])>0 else stdshift_penalty for a,b in zip(breaks_x[max(-turn,0):],breaks_y[max(turn,0):])]) for turn in  xrange(minisize-len(breaks_x),len(breaks_y)-minisize+1)]
		#differences_reverse=[sum([ calculatedistance(a[0],b[0]) + calculatedistance(a[1],b[1]) if (a[1]-a[0])*(b[1]-b[0])>0 else stdshift_penalty for a,b in zip(reverse_x[max(-turn,0):],breaks_y[max(turn,0):])]) for turn in  xrange(minisize-len(reverse_x),len(breaks_y)-minisize+1)]
		
		if min(differences_reverse)<min(differences):
			std=-1
			indent=differences_reverse.index(min(differences_reverse))+(minisize-len(breaks_x))
		else:
			std=1
			indent=differences.index(min(differences))+(minisize-len(breaks_x))
		"""
		
		differences=sum([ calculatedistance(a[0],b[0]) + calculatedistance(a[1],b[1]) if (a[1]-a[0])*(b[1]-b[0])>0 else calculatedistance(a[1],b[0]) + calculatedistance(a[0],b[1])+stdshift_penalty for a,b in zip(breaks_x,breaks_y)]) 
		
		differences_reverse=sum([ calculatedistance(a[0],b[0]) + calculatedistance(a[1],b[1]) if (a[1]-a[0])*(b[1]-b[0])>0 else calculatedistance(a[1],b[0]) + calculatedistance(a[0],b[1])+stdshift_penalty for a,b in zip(reverse_x,breaks_y) ]) 
		
		if differences_reverse<differences:
			std=-1
		else:
			std=1
	
		return min(differences,differences_reverse),indent,std


def find_alignments(x,y):
	
	
	front_x=list(x[0])
	back_x=list(x[-1])
	
	front_y=list(y[0])
	back_y=list(y[-1])
	
	if x[0][1]>=x[0][0]:
		front_x[0]=0
	else:
		front_x[0]=1000000
	
	if y[0][1]>=y[0][0]:
		front_y[0]=0
	else:
		front_y[0]=1000000
	
	if x[-1][1]>=x[-1][0]:
		back_x[1]=1000000
	else:
		back_x[1]=0
	
	if y[-1][1]>=y[-1][0]:
		back_y[1]=1000000
	else:
		back_y[1]=0
	
	breaks_x=[tuple(front_x)]+[a for a in x[1:-1] if a[1]>=0]+[tuple(back_x)]
	breaks_y=[tuple(front_y)]+[a for a in y[1:-1] if a[1]>=0]+[tuple(back_y)]

	x_findalignment=[]
	y_findalignment=[]
	
	for index,x_alignment in enumerate(x):
		
		
		similarity_f=[(calculatedistance(x_alignment[0],b[0]) + calculatedistance(x_alignment[1],b[1]) if (x_alignment[1]-x_alignment[0])*(b[1]-b[0])>0 else stdshift_penalty) if i not in y_findalignment else 100000 for i,b in enumerate(y)]
		
		x_reverse=makereverse([x_alignment])[0]
		similarity_r=[ (calculatedistance(x_reverse[0],b[0]) + calculatedistance(x_reverse[1],b[1]) if (x_reverse[1]-x_reverse[0])*(b[1]-b[0])>0 else stdshift_penalty) if i not in y_findalignment else 100000 for i,b in enumerate(y) ]
		
		if min(similarity_f+similarity_r)<100:
			
			if min(similarity_f)<=min(similarity_r):
				y_findalignment.append(similarity_f.index(min(similarity_f+similarity_r)))
			else:
				y_findalignment.append(similarity_r.index(min(similarity_f+similarity_r)))
			
			x_findalignment.append(index)
	
	alignments=[[] if i not in x_findalignment else y_findalignment[x_findalignment.index(i)] for i in xrange(len(x))]
	
	return alignments
			




def indent_alignments(cluster,names):
	
	
	sort_cluster_index=sorted(range(len(cluster)),key=lambda x:len(cluster[x]),reverse=True)
	
	sort_cluster=[cluster[x] for x in sort_cluster_index]
	sort_names=[names[x] for x in sort_cluster_index]
	
	longest=sort_cluster[0]
	
	indents=[0]
	stds=[1]
	for molecule in sort_cluster[1:]:
		
		
		diff,indent,std=similarity(molecule,longest,1)
		
		indents.append(indent)
		stds.append(std)
	
	indents=[x-min(indents) for x in indents]

	sort_cluster=[[indent]+molecule if std>0 else [indent]+makereverse(molecule) for indent,molecule,std in zip(indents,sort_cluster,stds)]
	
	
	return sort_cluster,sort_names


def find_cmapindex(IDs,left_end, right_end):
	
	target_region={rid:pos for rid,pos in IDs.items() if pos>=left_end and pos<=right_end}

	return min(target_region),max(target_region)



def read_cmap(qmapfile):
	
	cmap=pd.read_csv(qmapfile, sep='\s+', comment='#', names=readheader(qmapfile))
	contigIDs=map(int,list(cmap['CMapId']))
	Positions=list(cmap['Position'])
	
	allcontigs=cl.defaultdict(list)
	for contigid, pos in zip(contigIDs, Positions):
			
		allcontigs[contigid].append(pos)

	return dict(allcontigs)


def breakalignments(alignment,rmapIDtoposi,qmapIDtoposi):
	
	if not len(alignment):
		return alignment
		
	alignments_out=[]
	new_chunk=[]
	last_index=alignment[0][0]
	last_qindex=alignment[0][1]
	last_posi=rmapIDtoposi[last_index]
	last_qposi=qmapIDtoposi[last_qindex]
	
	for pair in alignment[1:]:
		
		current_index=pair[0]
		current_qindex=pair[1]
		current_posi=rmapIDtoposi[current_index]
		current_qposi=qmapIDtoposi[current_qindex]
		
		if (abs(current_posi-last_posi)>10000 and abs(current_index-last_index)>2) or (abs(current_qposi-last_qposi)>10000 and abs(current_qindex-last_qindex)>2):
			
			alignments_out.append(new_chunk)
			
			new_chunk=[pair]
		else:
			new_chunk.append(pair)
			
		last_posi=current_posi
		last_index=current_index
		last_qposi=current_qposi
		last_qindex=current_qindex
		
	alignments_out.append(new_chunk)
	
	return alignments_out	



def locate_qmap(alignment, strand, rmapIDtoposi,qmapIDtoposi,left_coordi,right_coordi,extend_size_up,extend_size_down):
	
	if strand=='+':
		strand=1
	else:
		strand=-1
	
	ref_left_cut=left_coordi-extend_size_up
	ref_right_cut=right_coordi+extend_size_down
		
	anchored_region=[]
	for pair in alignment:
		rcoordi=rmapIDtoposi[pair[0]]
		qcoordi=qmapIDtoposi[pair[1]]
		
		if rcoordi>=ref_left_cut and rcoordi<=ref_right_cut:
			anchored_region.append((rcoordi,qcoordi))
			
	anchored_region.sort()
	
	if not len(anchored_region):
		return -1,-1
	
	if extend_size_down:
		
		extend_size=ref_right_cut-anchored_region[-1][0]
	
		query_left=anchored_region[0][1]
		query_right=anchored_region[-1][1]+strand*extend_size
	
	if extend_size_up:
		
		extend_size=anchored_region[0][0]-ref_left_cut
		
		query_left=anchored_region[0][1]-strand*extend_size
		query_right=anchored_region[-1][1]	
	
	qcut_left=min(query_left,query_right)
	qcut_right=max(query_left,query_right)
	
	qstart=-1
	for i, pos in enumerate(qmapIDtoposi):
		
		if pos>=qcut_left:
			qstart=i
			break
	
	qend=0
	for i, pos in enumerate(qmapIDtoposi[qstart:]):
			
		qend=i
		if pos>qcut_right:
			break

	qend+=qstart

	return qstart, qend
	
def subset_byanchor(alignment, query_cuts):
	
	query_left,query_right=query_cuts
	
	if query_left<0 or not len(alignment):
		return []
		
	result=[pair for pair in alignment if pair[1]>=query_left and pair[1]<=query_right]
		

	return result



def readxmap(xmapfile,rmapIDtoposi,qmapIDtoposi,left_coordi,right_coordi,extend_size_up,extend_size_down):
		
	xmap=pd.read_csv(xmapfile, sep='\s+', comment='#', names=readheader(xmapfile))
	
	alignments=list(xmap['Alignment'])
	
	allcontigs=list(xmap['QryContigID'])
	
	allstrands=list(xmap['Orientation'])
	
	target_alignments=[[tuple([int(a) for a in x.split(',')]) for x in re.findall(r'(?<=\()(\d+\,\d+)(?=\))',alignment)] for alignment in alignments]
	
		
	out_alignments, out_allcontigs, out_allstrands=[],[],[]
	
	contigs=cl.defaultdict(list)
	contigs_target={}
	
	query_cuts={}
	for alignment, contig, strand in zip(target_alignments, allcontigs, allstrands):

		new_qstart, new_qend=locate_qmap(alignment, strand, rmapIDtoposi,qmapIDtoposi[contig],left_coordi,right_coordi,extend_size_up,extend_size_down)

		qstart,qend=query_cuts.get(contig,(-1,-1))
		
		if new_qstart>=0:
			
			if qstart<0:
				query_cuts[contig]=(new_qstart,new_qend)
			else:
				query_cuts[contig]=(min(qstart,new_qstart),max(new_qend, qend))

	for alignment, contig, strand in zip(target_alignments, allcontigs, allstrands):
		
		alignment =subset_byanchor(alignment,query_cuts.get(contig,(-1,-1)))
		
		if strand=='+':
			alignments_cuts=[x for x in breakalignments(alignment,rmapIDtoposi,qmapIDtoposi[contig]) if len(x)>2]
		else:
			alignments_cuts=[x for x in breakalignments(alignment[::-1],rmapIDtoposi,qmapIDtoposi[contig]) if len(x)>2]
		
		contigs[contig].extend(alignments_cuts)
		
	contigs={contig:sorted(value, key=lambda x:min(x[0][1],x[-1][1])) if len(value) else [] for contig, value in contigs.items()}#sorted based query
	
	new_contigs={}
	for name, alignments in contigs.items():
	
		contig,starts=summary_contig(rmapIDtoposi,qmapIDtoposi[name],[(x[0][1],x[-1][1],x[0][0],x[-1][0]) for x in alignments], left_coordi,right_coordi)
		
		if len(contig)==1:
			new_contigs[str(name)]=contig[0]
			
		else:
			
			for i,contig0 in enumerate(contig):
				
				new_contigs[str(name)+"-"+str(starts[i])]=contig0
	
	return new_contigs
	
def summary_contig(rmapIDtoposi,qmapIDtoposi, breaks,left_coordi,right_coordi):
	
	if not len(breaks):
		return [],[]
	
	alignment_start=min(rmapIDtoposi[breaks[0][3]], rmapIDtoposi[breaks[-1][3]])
	alignment_end=max(rmapIDtoposi[breaks[0][3]], rmapIDtoposi[breaks[-1][3]])
	
	
	if max(alignment_end, right_coordi)-min(alignment_start, left_coordi)-(alignment_end-alignment_start)-(right_coordi-left_coordi)>0:
		return [],[]
	
	splits=[]
	new_breaks=[breaks[0][2:]]
	starts=[qmapIDtoposi[breaks[0][0]]]
	
	current_index=-1
	old_end=breaks[0][1]
	for current_index,break0 in enumerate(breaks[1:]):
		
		q_start,q_end,r_start,r_end=break0
		
		if old_end-q_end>10000:
			
			splits.append(new_breaks)
			new_breaks=[(r_start,r_end)]
			starts.append(qmapIDtoposi[q_start])
			continue
		
		if q_start-old_end>1:
			new_breaks.append((-1,int(abs(qmapIDtoposi[q_start]-qmapIDtoposi[old_end]))))
			
		new_breaks.append((r_start,r_end))
		
		old_end=q_end

	splits.append(new_breaks)
	
	
	return splits, starts


def readheader(filename):
	
	thefile=open(filename, mode='r')
	
	header=''
	for line in thefile:
		
		split=line.split()
		if split[0]=="#h":
			
			header=split[1:]
			break
			
	thefile.close()
	
	return header
	

class cluster:
	
	def __init__(self):
		
		self.clusters={}
	
	def load_data(self,xmapfile, qmapfile, rmapfile):
		
		self.xmapfile, self.qmapfile, self.rmapfile=xmapfile, qmapfile, rmapfile
	
	def subset(self,chr0, anchorstart, anchorend, anchor_extend_up, anchor_extend_down):
		
		self.chr=int(chr0)
		self.start=anchorstart
		self.end=anchorend
		self.anchor_extend_up=anchor_extend_up
		self.anchor_extend_down=anchor_extend_down
		
		ref_cmap=pd.read_csv(self.rmapfile, sep='\s+', comment='#', names=readheader(self.rmapfile))
		ref_cmap=ref_cmap[(ref_cmap['CMapId']==self.chr)]
		
		SiteIDs=list(ref_cmap['SiteID'])
		Positions=list(ref_cmap['Position'])

		self.IDtoposi={id0:int(posi) for id0,posi in zip(SiteIDs,Positions)}
				
		contigs= readxmap(self.xmapfile,self.IDtoposi,read_cmap(self.qmapfile),self.start, self.end, self.anchor_extend_up, self.anchor_extend_down)

		self.contig=contigs
		
	def cluster(self):
			
		self.clusters={}
		for contig,breaks in self.contig.items():
			
			if not len(breaks):
				continue
			
			intersections=[]
			for names,cluster in self.clusters.items():
				
				names=names.split('_')
				intersections.append({names[i]:1.0*similarity(breaks, x) for i,x in enumerate(cluster)})
				


			if len(intersections)>0:
				
				pass_filter_index=[i for i,x in enumerate(intersections) if 1.00*sum([y**2 for y in x.values()])/len(x.values()) <cluster_max]
				
			
			else:
				pass_filter_index=[]

				
			
			if len(pass_filter_index)>0:
				
				allindex_names=[self.clusters.keys()[i] for i in pass_filter_index]
				
				sample_name=str(contig).split('000')[0].split('-')[0]
						
				best_cluster_index=sorted(pass_filter_index,key=lambda x: -10000*self.clusters.keys()[x].count(sample_name)+1.00*sum(intersections[x].values())/len(intersections[x].values())-len(intersections[x].values())*0.1)[0]
						
				index_name=self.clusters.keys()[best_cluster_index]
				
				index_name_new=index_name+'_'+str(contig)
				
				self.clusters[index_name_new]=self.clusters[index_name]+[breaks]
				
				del self.clusters[index_name]
				
			
			else:
				
				self.clusters[str(contig)]=[breaks]
	
		self.clusters={name:data for name,data in self.clusters.items() if len(data)>1}
		
	def output(self, output):
		
		gap=''.join([' ' for x in xrange(12)])
		gap2=''.join([' ' for x in xrange(22)])

		self.out=[]
		out2=[]
		
		spaces=[''.join([" " for x in xrange(l)]) for l in xrange(100)]
				
		for index, (name, cluster) in  enumerate(self.clusters.items()):
								
			names=name.split('_')
			
			cluster_sort_key=sorted(range(len(names)), key=lambda x: names[x])
			cluster=[cluster[i] for i in cluster_sort_key]
			names=sorted(names)
						
			cluster_indented,names=indent_alignments(cluster, names)
			
				
			out0='\ncluster_%d\n'%index
			out02='\ncluster_%d\n'%index
			for molecule,name in zip(cluster_indented,names):
				
				indent=molecule[0]
				molecule=molecule[1:]
				
				text_alignments=['Insert:{:d}'.format(x[1]) if x[0]<0 else '({:d}-{:d})+'.format(x[0],x[1]) if x[1]>=x[0] else '({:d}-{:d})-'.format(x[0],x[1]) for x in molecule]
				
				text_alignments=[spaces[(22-len(x))/2]+x+spaces[(23-len(x))/2] for x in text_alignments]
				
				out0=out0+name+','+','.join([gap for x in xrange(indent)]+text_alignments)+' \n'	
				
				text_alignments2=['Insert:{:d}'.format(x[1]) if x[0]<0 else '({:d}-{:d})+'.format(self.IDtoposi[x[0]],self.IDtoposi[x[1]]) if x[1]>=x[0] else '({:d}-{:d})-'.format(self.IDtoposi[x[0]],self.IDtoposi[x[1]]) for x in molecule]
				text_alignments2=[spaces[(22-len(x))/2]+x+spaces[(23-len(x))/2] for x in text_alignments2]
				
				out02=out02+name+spaces[12-len(name)]+','+','.join([gap2 for x in xrange(indent)]+text_alignments2)+' \n'	
				
				
			self.out.append(out0)
			out2.append(out02)
			
		self.out='\n'.join(sorted(self.out,key=lambda x: len(x), reverse=True))
		
		out2='\n'.join(sorted(out2,key=lambda x: len(x), reverse=True))
		with open(output,mode='w') as f:
			f.write(out2)
		f.close()
		
	
	def summary(self, smfile):
		
		cluster_references=[]

		cluster_names=re.findall(r'cluster_\d+',out)
		all_clusters=self.out.split('cluster')[1:]


		for cluster in all_clusters:

			alllines=zip(*[x.split(',')[1:] for x in cluster.splitlines() if len(x.strip())>0 and ',' in x])
			allpieces=[list([[int(a) for a in re.findall(r'\d+',y)] for y in x]) for x in alllines]
			
			reference=[]
			for piece in allpieces:

				front=[x[0] for x in piece if len(x)>1]
				end=[x[1] for x in piece if len(x)>1]
				
				front_common=sorted(list(set(front)),key=lambda x:front.count(x), reverse=True)[0]
				end_common=sorted(list(set(end)),key=lambda x:front.count(x), reverse=True)[0]
				reference.append((front_common,end_common))

			cluster_references.append(reference)


		longest=sorted(cluster_references,key=lambda x:len(x))[-1]

		out=[]
		for name,reference in zip(cluster_names,cluster_references):
			
			
			alignment=find_alignments(reference,longest)

			out0=name+':\t'+' , '.join([gap if len(x)<2 else '({:d}-{:d})+'.format(x[0],x[1]) if x[1]>=x[0] else '({:d}-{:d})-'.format(x[0],x[1]) for x in reference])

			out.append(out0)
			
		out='\n'.join(out)

		cluster_names,cluster_references_indented=indent_alignments(cluster_references,cluster_names)
		out='\n'.join([name+':\t'+' , '.join(['({:d}-{:d})+'.format(self.IDtoposi[x[0]],self.IDtoposi[x[1]]) if x[1]>=x[0] else '({:d}-{:d})-'.format(self.IDtoposi[x[0]],self.IDtoposi[x[1]]) for x in reference]) for name,reference in zip(cluster_names,cluster_references)])

		with open(smfile, mode='w') as f:
			f.write(out)
		f.close()
			


def run(parser):
	
	global stdshift_penalty, cluster_max, base_penalty
	
	stdshift_penalty, cluster_max, base_penalty=parser.stdshift_penalty, parser.cluster_max, parser.base_penalty
	
	thecluster=cluster()

	thecluster.load_data(parser.xmap,parser.qmap, parser.rmap)
	
	thecluster.subset(parser.chr,parser.start,parser.end,parser.upstream_extend,parser.downstream_extend)
	
	thecluster.cluster()

	thecluster.output(parser.out)
		

def main():
	parser=argparse.ArgumentParser(description="script to cluster bionano local alignments")
	parser.add_argument("-x","--xmap",help="xmap file" ,dest="xmap",type=str)
	parser.add_argument("-q","--query_cmap",help="query_cmap file" ,dest="qmap",type=str)
	parser.add_argument("-r","--ref_cmap",help="rmap file",dest="rmap",type=str)
	parser.add_argument("-o","--output",help="output file",dest="out",type=str, default="./output")
	parser.add_argument("-c","--chr",help="target region chromosome" ,dest="chr",type=int)
	parser.add_argument("-s","--start",help="target region start" ,dest="start",type=int)
	parser.add_argument("-e","--end",help="target region end" ,dest="end",type=int)
	parser.add_argument("-u","--upstream_extend",help="extending anchor region upstream size" ,dest="upstream_extend",type=float, default=0)
	parser.add_argument("-d","--downstream_extend",help="extending anchor region downstream size" ,dest="downstream_extend",type=float, default=0)
	parser.add_argument("-m","--cluster_max",help="maximum difference to cluster" ,dest="cluster_max",type=float, default=4000)
	parser.add_argument("-p","--stdshift_penalty",help="strandshift_penalty" ,dest="stdshift_penalty",type=float, default=4000)
	parser.add_argument("-b","--base_penalty",help="defaulty baseline penalty" ,dest="base_penalty",type=float, default=4000)
	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)	
	

if __name__=='__main__':
	
	main()



