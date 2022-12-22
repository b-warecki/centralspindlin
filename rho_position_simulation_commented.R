#associated with Figure S1 for "Cleavage plane positioning by centralspindlin-mediated RhoGEF transport" (Warecki and Tao, 2023)
#please direct questions to: bwarecki@ucsc.edu
#this simulation will determine the distribution of Centralspindin-RhoGEF molecules along a microtubule over time

#load required packages
library(ggplot2)
library(dplyr)
library(RColorBrewer)

#simulation function:
RhoGEF.pos.time <- function(total.time, numCS, Plk.time, retention, relative.velocity) {
		#total.time=length of time for simulation to run
		#numCS=number of CS molecules to add to simulation
		#Plk.time=time after the start of simulation at which point RhoGEF can bind CS (represents polo phosphorylation)
		#retention=the retention rate of cargo-less CS at plus-end tips prior to polo phosphorylation
		#relative.velocity=how much faster/slower CS-RhoGEF transport is than cargo-less CS walking along microtubule
	time <- seq(0, (total.time-1)) #creates a sequence of time points
	MT.length <- 60 #arbitrary length of microtubule (position 60=plus-end tip)
		#below to name output variables
	rhogef.position <- vector() 
	rhogef.time <- vector()
	CS.num <- vector()
	CS.number <- 1 #identifier for CS molecule
		#below to define dataframes for simulation output
	d.f <- vector()
	for (i in 1:numCS) { #for every CS molecule
		CS.position <- sample(seq(1,MT.length), 1, replace=TRUE) #CS randomly associates along length of microtubule
			#below defines starting values that will be changed throughout simulation
		RhoGEF.CS.new <- 0
		RhoGEF.CS <- 0
		RhoGEF.t <- 0
		for (i in seq_along(time)) { #for every timepoint
			if (i < Plk.time) { #before polo phosphorylation
				CS.position <- CS.position+1 #CS advances forward 1 position
				if (CS.position>=MT.length) { #if CS reaches plus-end tip
					CS.ret <- sample(0:1, prob=c((1-retention), retention))[1] #determines if cargo-less CS remains on the plus-end tip based on entered retention value
					if (CS.ret==1) {
						CS.position <- MT.length #cargo-less CS is retained on plus-end tip
					} else {
						CS.position <- sample(seq(0,MT.length), 1, replace=TRUE) #CS dissociates and rebinds randomly along length of microtubule
					}
				}
				RhoGEF.t <- RhoGEF.t+1 #advances time counter for output
				RhoGEF.pos <- 0 #no RhoGEF on microtubule prior to Polo phosphorylation
				rhogef.time <- append(rhogef.time, RhoGEF.t) #add current timepoint to output
				rhogef.position <- append(rhogef.position, RhoGEF.pos) #add current RhoGEF position to output
				CS.num <- append(CS.num, CS.number) #add current CS molecule to output
			} else { #after polo phosphorylation
				if (RhoGEF.CS==0) { #for any cargo-less CS molecule
					RhoGEF.CS.new <- sample(0:1, prob=c(0.1,0.9))[1] #determines if RhoGEF binds CS		
				}
				if (RhoGEF.CS.new==1) {
					RhoGEF.CS <- 1 #newly bound CS will be treated as a CS-RhoGEF complex
				}
				if (RhoGEF.CS==1) {
					CS.position <- CS.position+(1*relative.velocity) #advances bound CS microtubule based on entered relative.velocity value
				} else {
					CS.position <- CS.position+1 #cargo-less CS advances forward 1 position
				}
				if (CS.position>=MT.length && RhoGEF.CS==1) { #for CS-RhoGEF complex at plus-end tip
					CS.ret <- sample(0:1, prob=c(0.08,0.92))[1] #probability CS-RhoGEF is retained at plus-end tip
					if (CS.ret==1) {
						CS.position <- MT.length #CS-RhoGEF retained at plus-end tip
					} else {
						CS.position <- sample(seq(0,MT.length), 1, replace=TRUE) #CS-RhoGEF dissociates and binds randomly along length of microtubule
					}
				}
				if (CS.position>=MT.length && RhoGEF.CS==0) { #for cargo-less CS at plus-end tip
					CS.ret <- sample(0:1, prob=c((1-retention),retention))[1] #probability CS is retained at plus-end tip based on entered retention value
					if (CS.ret==1) {
						CS.position <- MT.length #cargo-less CS retained at plus-end tip
					} else {
						CS.position <- sample(seq(0,MT.length), 1, replace=TRUE) #CS dissociates and binds randomly along length of microtubule
					}
				}
				RhoGEF.t <- RhoGEF.t+1 #advances time counter for output
				if (RhoGEF.CS==1) { #for CS-RhoGEF complexes
					RhoGEF.pos <- CS.position #position of Rho on microtubule = position of CS on microtubule
				} else {
					RhoGEF.pos <- 0 #for cargo-less CS, no RhoGEF position on microtubule
				}
				rhogef.time <- append(rhogef.time, RhoGEF.t) #add current timepoint to output
				rhogef.position <- append(rhogef.position, RhoGEF.pos) #add current Rho position to output
				CS.num <- append(CS.num, CS.number) #add current CS molecule to output
			}
			d.f <- data.frame(Pos=rhogef.position, Time=rhogef.time, CS=CS.num)	#creates data frame with each RhoGEF position along the microtubule at each time for each CS molecule
		}
		CS.number <- CS.number + 1 #create new identifier for next CS molecule
	}
	d.f.afterpolo <- filter(d.f, Time>=25) #filter data for RhoGEF position after polo phosphorylation
	return(d.f.afterpolo)
}

#run simulations under different conditions
test.null <- RhoGEF.pos.time(50,1000,25,0.0000001,1) #no cargo-less retention + no CS stimulation by RhoGEF
test.retain <- RhoGEF.pos.time(50,1000,25,0.89,1) #cargo-less retention + no CS stimulation by RhoGEF
test.retain_transport <- RhoGEF.pos.time(50,1000,25,0.89,(11/6)) #cargo-less retention + CS stimulation by RhoGEF

#for graphical output
#create histograms to place CS-RhoGEF molecules in bins of 5 positions along microtubule
histogram.rearrange.5 <- function(data) {
	max.time <- max(data$Time) #determine maximum time value in simulation
	initial.time <- 24
	data.binned.total <- vector() #define output
	for (i in 24:max.time) {
		d.f <- filter(data, Time==initial.time) #pull out all positions for a given time
		hist.count <- hist(d.f$Pos, plot=FALSE, breaks=seq(0,60,by=5))$counts #group counts into bins of 5
		hist.break <- seq(0,60,by=5)[2:13] #label bins
		time.list <- rep(initial.time, 12) #label time
		data.binned <- data.frame(count=hist.count, pos=hist.break, time.l=time.list) #create data frame for each time point
		data.binned.total <- rbind(data.binned.total, data.binned) #combine data frames from all time points
		initial.time <- initial.time+1 #advance time counter
	}
	return(data.binned.total)
}

#create binned data for simulations
test.null.rearrange <- histogram.rearrange.5(test.null)
test.retain.rearrange <- histogram.rearrange.5(test.retain)
test.retain_transport.rearrange <- histogram.rearrange.5(test.retain_transport)

#qualitative distribution of RhoGEF along a microtubule over time (Figure S1C)
dev.new(width=2, height=2) 
ggplot(test.null.rearrange, aes(x=pos, y=time.l, fill=count)) +
	geom_tile() +
	xlim(6,65) +
	ylim(24,41) +
	scale_fill_distiller(palette="Spectral", limits=c(0,600)) +
	theme_classic() +
	annotate("segment", x=6, xend=65, y=24.5, yend=24.5, linetype="dashed") +
	guides(fill=FALSE) +
	labs(y="Time (abu)",
		 x="Position along MT") +
	theme(axis.text.x=element_blank())
ggsave("rhogef.null.tile.jpeg", dpi=1000)
	
dev.new(width=2, height=2) 
ggplot(test.retain.rearrange, aes(x=pos, y=time.l, fill=count)) +
	geom_tile() +
	xlim(6,65) +
	ylim(24,41) +
	scale_fill_distiller(palette="Spectral", limits=c(0,600)) +
	theme_classic() +
	annotate("segment", x=6, xend=65, y=24.5, yend=24.5, linetype="dashed") +
	guides(fill=FALSE) +
	labs(y="Time (abu)",
		 x="Position along MT") +
	theme(axis.text.x=element_blank())
ggsave("rhogef.retain.tile.jpeg", dpi=1000)
	
dev.new(width=2, height=2) 
ggplot(test.retain_transport.rearrange, aes(x=pos, y=time.l, fill=count)) +
	geom_tile() +
	xlim(6,65) +
	ylim(24,41) +
	scale_fill_distiller(palette="Spectral", limits=c(0,600)) +
	theme_classic() +
	annotate("segment", x=6, xend=65, y=24.5, yend=24.5, linetype="dashed") +
	guides(fill=FALSE) +
	labs(y="Time (abu)",
		 x="Position along MT") +
	theme(axis.text.x=element_blank())
ggsave("rhogef.retain_transport.jpeg", dpi=1000)
	
#quantitative comparison of RhoGEF at microtubule plus-end (last 5 positions of microtubule) for different simulations (Figure S1D)
test.null.rearrange.60 <- filter(test.null.rearrange, pos==60) #pull out data for microtubule end
test.null.rearrange.60 <- mutate(test.null.rearrange.60, test="null") #define test
test.retain.rearrange.60 <- filter(test.retain.rearrange, pos==60) #pull out data for microtubule end
test.retain.rearrange.60 <- mutate(test.retain.rearrange.60, test="retain") #define test
test.retain_transport.rearrange.60 <- filter(test.retain_transport.rearrange, pos==60) #pull out data for microtubule end
test.retain_transport.rearrange.60 <- mutate(test.retain_transport.rearrange.60, test="retain_transport") #define test
test.60 <- rbind(test.null.rearrange.60, test.retain.rearrange.60, test.retain_transport.rearrange.60) #combine dataframes

dev.new(width=2.25, height=2.25)
ggplot(test.60, aes(x=time.l, y=count, fill=rev(test))) +
	geom_line(size=1.5) +
	geom_area(position="identity", alpha=0.3) +
	theme_classic() +
	xlim(24,40) +
	guides(fill=FALSE) +
	theme(axis.text.x=element_text(size=11),
		  axis.text.y=element_text(size=11)) +
	scale_fill_manual(values=c("grey", "blue", "magenta")) +
	labs(x="Time (abu)",
		 y="RhoGEF# at plus-end") +
	annotate("segment", x=24, xend=24, y=0, yend=Inf, linetype="dashed")
ggsave("rhogef.atplustend.jpeg", dpi=1000)
	
#quantitative comparison of RhoGEF distribution along the microtubule at time=40 for each simulation (Figure S1E)
test.null.40 <- filter(test.null, Time==40) #pull RhoGEFF positions for time=40
test.null.40 <- mutate(test.null.40, test="null") #define test
test.retain.40 <- filter(test.retain, Time==40) #pull RhoGEFF positions for time=40
test.retain.40 <- mutate(test.retain.40, test="retain") #define test
test.retain_transport.40 <- filter(test.retain_transport, Time==40) #pull RhoGEFF positions for time=40
test.retain_transport.40 <- mutate(test.retain_transport.40, test="retain_transport") #define test
test.40 <- rbind(test.null.40, test.retain.40, test.retain_transport.40) #combine data frames

dev.new(width=2, height=2.25)
ggplot(test.40, aes(x=Pos, group=test, fill=test)) +
	geom_density(size=1, alpha=0.5) +
	theme_classic() +
	guides(fill=FALSE) +
	scale_fill_manual(values=c("magenta", "blue", "grey")) +
	labs(x="Position along MT",
		 y="Density of RhoGEF\nat Time=40")
ggsave("rhogef.dist.time40.jpeg", dpi=1000)
