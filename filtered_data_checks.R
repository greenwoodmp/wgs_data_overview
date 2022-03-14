rm(list=ls())

#------------------------#
#### Loading Packages ####
#------------------------#

#install the pacman package manager if necessary
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")

pacman::p_load(here, # to provide document paths relative to the R project
               ggplot2, # for plotting
               dplyr, # for easy dataframe reconstruction
               reshape,# for the melt function
               tidyr,
               ggpubr) 

#-------------------------------------------------#
#### Assessing Reads Following Trimmomatic Run ####
#-------------------------------------------------#

dat <- read.csv(here("./Trimmed Read Statistics/trimming_stats.csv"), header=FALSE)

# The names are given based on the order in which data was assembled in the 1.3-TrimmingStats.sh bash script
names(dat) <- c("Accession", "Raw", "Paired", "Unpaired_Forward", "Unpaired_Reverse")

dat$Forward_Lost <- dat$Raw - dat$Paired - dat$Unpaired_Forward
dat$Reverse_Lost <- dat$Raw - dat$Paired - dat$Unpaired_Reverse
dat$Total_Lost <- dat$Forward_Lost + dat$Reverse_Lost

# Sweep will divide entries ina column by a vector. Here, we use the parameter
#   "1" for each row to be divided by the raw number of reads in that row
percent_dat <- sweep(dat[,-1], 1, as.numeric(dat$Raw),"/")*100
percent_dat <- data.frame(dat$Accession, percent_dat)
names(percent_dat)[1] <- c("Accession")
levels(percent_dat) <- as.factor(colnames(percent_dat))

#-------------------------------------#
#### Boxplots of the Retained Data ####
#-------------------------------------#

# Melt the dataframe for plotting purposes
melted <- gather(percent_dat[,-2],# having raw reads is pointless 
                 "Readset",
                 "Reads",
                 -Accession)
melted$Readset <- factor(melted$Readset,
                         levels=c("Accession", "Raw", "Paired", "Unpaired_Forward", "Unpaired_Reverse","Forward_Lost", "Reverse_Lost", "Total_Lost"))

means <- aggregate(Reads ~ Readset, melted, mean)
means$Reads <- round(means$Reads, 2)

p <- ggplot(melted, aes(x=as.factor(Readset), y=Reads)) +
  geom_boxplot()+
  geom_text(data = means, aes(label = Reads, y = Reads + 5))+
  theme_bw()+
  ylab("Reads as % of Total Raw Reads Per Accession")+
  xlab("Trimmed Read Groups")


#-------------------------------------------------#
#### Violin Plots of Estimated Genome Coverage ####
#-------------------------------------------------# 

# Estimate the total genome coverage of reads based on the assembled BST1 genome
#   size of 266544826 bp (obtained on the cluster using "cat /bettik/PROJECTS/pr-mosquito/COMMON/Coenonympha/Genomes/BST1.hifiasm.p_ctg.cleaned.sorted.fa | paste - - | cut -f2 | tr -d '\n' | wc -c")
#   and an average read size of 150bp
paired_coverage <- (dat$Paired*2*150)/533088877 # dat paired is multipied by 2 to include reverse and forward reads
total_coverage <- ((dat$Paired*2+dat$Unpaired_Forward+dat$Unpaired_Reverse)*150)/533088877 # dat paired is multipied by 2 to include reverse and forward reads
best_coverage <- (dat$Raw*2*150)/533088877 # dat Raw_Forward is multipied by 2 to include reverse and forward reads; this is a scenario in which no reads were lost

cov <- data.frame(dat$Accession, paired_coverage, total_coverage, best_coverage)
names(cov)[1] <- "Accession"
meltcov <- gather(cov,# having raw reads is pointless 
                 "Readset",
                 "Theoretical_Coverage",
                 -Accession)
meltcov$Readset <- as.factor(meltcov$Readset) 

covmeans <- aggregate(Theoretical_Coverage ~ Readset, meltcov, median)

dp <- ggplot(meltcov, aes(x=Theoretical_Coverage, color=Readset, fill=Readset))+
  geom_density()+
  geom_vline(data=covmeans, aes(xintercept=Theoretical_Coverage, color=Readset),
             linetype="dashed")+
  facet_grid(. ~ Readset)

paired_mean <- data.frame(mean(cov$paired_coverage))
paired_sd <- data.frame(sd(cov$paired_coverage))

# Make density curve for plotting colours/points under the curve (see https://stackoverflow.com/questions/20355849/ggplot2-shade-area-under-density-curve-by-group)
dc <- data.frame(list(x=density(cov$paired_coverage)$x, y=density(cov$paired_coverage)$y))

# Get quantiles of the data which we consider to be at the extreme ends of the distribution
lowerq <- quantile(cov$paired_coverage,0.25)
upperq <- quantile(cov$paired_coverage,0.25)

dens <- ggplot(cov, aes(x=paired_coverage))+
  geom_density(size=1, fill="#87BBA2")+
  geom_vline(data=paired_mean, aes(xintercept=mean.cov.paired_coverage.),
             linetype="dashed", size=1.5, colour='#1E3F20')+
  # Use ribbon and the dc distribution to highlight the curve outside 1 SD of the data
  geom_ribbon(data=subset(dc,x<(paired_mean$mean.cov.paired_coverage.-paired_sd$sd.cov.paired_coverage.)),
              aes(x=x,ymax=y),ymin=0,fill="#D34E24", alpha=0.8)+
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background=element_rect(colour='black', fill='white')
  )+
  xlab("Theoretical Genome Coverage")+
  ylab("Density")+
  geom_text(data=paired_mean, aes(x=mean.cov.paired_coverage., y=0.2, label=round(mean.cov.paired_coverage.,2)))+
  geom_density(size=1) #This line is added to reestablish the density curve line over the orange part of the plot

### Actual Mapping Data

bwa <- read.csv(here("./Mapped Read Statistics/bwa_mapping_stats.csv"), header=FALSE)
names(bwa) <- c("Accession","Percent_Mapped","Percent_Mapped_Paired","Alignments","MAPQ>0","MAPQ>9","MAPQ>19","MAPQ>29","Mean_Coverage","Mean_Relative_Coverage","Breadth_Coverage")

dcr <- data.frame(list(x=density(bwa$Mean_Coverage)$x, y=density(bwa$Mean_Coverage)$y))
real_mean <- data.frame(mean(bwa$Mean_Coverage))
real_sd <- data.frame(sd(bwa$Mean_Coverage))

densreal <- ggplot(bwa, aes(x=Mean_Coverage))+
  geom_density(size=1, fill="#87BBA2")+
  geom_vline(aes(xintercept=mean(bwa$Mean_Coverage)),
             linetype="dashed", size=1.5, colour='#1E3F20')+
  # Use ribbon and the dc distribution to highlight the curve outside 1 SD of the data
  geom_ribbon(data=subset(dcr,x<(real_mean$mean.bwa.Mean_Coverage.-real_sd$sd.bwa.Mean_Coverage.)),
              aes(x=x,ymax=y),ymin=0,fill="#D34E24", alpha=0.8)+
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background=element_rect(colour='black', fill='white')
  )+
  xlab("BWA Genome Mapping Coverage")+
  ylab("Density")+
  geom_text(data=paired_mean, aes(x=mean(bwa$Mean_Coverage), y=0.25, label=round(mean(bwa$Mean_Coverage),2)))+
  geom_density(size=1)+ #This line is added to reestablish the density curve line over the orange part of the plot
  geom_point()

# Compare empirical total genome coverage to the theoretical genome coverage
together <- ggplot(bwa, aes(x=Mean_Coverage))+
  geom_histogram(data=cov, aes(x=paired_coverage, y=stat(density)), size=1, fill="#D34E24", alpha=0.8)+
  geom_density(data=cov, aes(x=paired_coverage), size=2, colour="#CC3000")+
  geom_vline(data=paired_mean, aes(xintercept=mean.cov.paired_coverage.),
             linetype="dashed", size=1.5, colour='#CC3000')+
  geom_text(data=paired_mean, aes(x=mean.cov.paired_coverage.+1.5, y=0.3, label=round(mean.cov.paired_coverage.,2)))+
  geom_histogram(aes(y=stat(density)), size=1, fill="#87BBA2", alpha=0.8)+
  geom_density(size=2, colour="#1E3F20")+
  geom_vline(aes(xintercept=mean(bwa$Mean_Coverage)),
             linetype="dashed", size=1.5, colour='#1E3F20')+
  # Use ribbon and the dc distribution to highlight the curve outside 1 SD of the data
  #geom_ribbon(data=subset(dcr,x<(real_mean$mean.bwa.Mean_Coverage.-real_sd$sd.bwa.Mean_Coverage.)),
              #aes(x=x,ymax=y),ymin=0,fill="#D34E24", alpha=0.8)+
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background=element_rect(colour='black', fill='white')
  )+
  xlab("BWA Genome Mapping Coverage")+
  ylab("Density")+
  geom_text(data=paired_mean, aes(x=mean(bwa$Mean_Coverage)+1.5, y=0.3, label=round(mean(bwa$Mean_Coverage),2)))

# Look at the distribution of genome coverage ONLY WITHIN regions that have mapped successfully
relcov <- ggplot(bwa, aes(x=Mean_Relative_Coverage))+
  geom_histogram(aes(y=stat(density)), size=1, fill="#87BBA2", alpha=0.8)+
  geom_density(size=2, colour="#1E3F20")+
  geom_vline(aes(xintercept=mean(bwa$Mean_Relative_Coverage)),
             linetype="dashed", size=1.5, colour='#1E3F20')+
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background=element_rect(colour='black', fill='white')
  )+
  xlab("BWA Genome Mapping Coverage in Successfully Mapped Regions")+
  ylab("Density")+
  geom_text(data=bwa, aes(x=mean(bwa$Mean_Relative_Coverage)+1.5, y=0.3, label=round(mean(bwa$Mean_Relative_Coverage),2)))

# Coverage the coverage across the genome to within regions that have mapped successfully
ggplot(bwa, aes(x=Mean_Coverage))+
  geom_histogram(aes(x=Mean_Relative_Coverage, y=stat(density)), size=1, fill="#5BB5EA", alpha=0.8)+
  geom_density(aes(x=Mean_Relative_Coverage), size=2, colour="#38799E")+
  geom_vline(aes(xintercept=mean(bwa$Mean_Relative_Coverage)),
             linetype="dashed", size=1.5, colour='#38799E')+
  geom_text(aes(x=mean(bwa$Mean_Relative_Coverage)+1.5, y=0.3, label=round(mean(bwa$Mean_Relative_Coverage),2)))+
  geom_histogram(aes(y=stat(density)), size=1, fill="#87BBA2", alpha=0.8)+
  geom_density(size=2, colour="#1E3F20")+
  geom_vline(aes(xintercept=mean(bwa$Mean_Coverage)),
             linetype="dashed", size=1.5, colour='#1E3F20')+
  # Use ribbon and the dc distribution to highlight the curve outside 1 SD of the data
  #geom_ribbon(data=subset(dcr,x<(real_mean$mean.bwa.Mean_Coverage.-real_sd$sd.bwa.Mean_Coverage.)),
  #aes(x=x,ymax=y),ymin=0,fill="#D34E24", alpha=0.8)+
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background=element_rect(colour='black', fill='white')
  )+
  xlab("BWA Genome Mapping Coverage")+
  ylab("Density")+
  geom_text(aes(x=mean(bwa$Mean_Coverage)+1.5, y=0.3, label=round(mean(bwa$Mean_Coverage),2)))+
  geom_segment(x=mean(bwa$Mean_Coverage)+0.5, xend=mean(bwa$Mean_Relative_Coverage)-0.5,
               y=0.25,yend=0.25,
               arrow = arrow(length = unit(0.5, "cm")),
               size=1.5,
               alpha=0.5,
               col="azure4")+
  geom_text(aes(x=((mean(bwa$Mean_Coverage)+mean(bwa$Mean_Relative_Coverage))/2),
                y=0.27),
            label="Dropping unmapped regions",
            hjust=0.5,
            size=3,
            col="azure4")



cor <- ggscatter(data=bwa, x="Mean_Coverage", y="Breadth_Coverage",
          add="reg.line",
          conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab="Mean Coverage",
          ylab="Portion of Genome Covered")

cor2 <- ggscatter(data=bwa, x="Mean_Relative_Coverage", y="Breadth_Coverage",
                 add="reg.line",
                 conf.int = TRUE,
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab="Mean Relative Coverage",
                 ylab="Portion of Genome Covered")

cor3 <- ggscatter(data=bwa, x="MAPQ>19", y="Mean_Relative_Coverage",
                 add="reg.line",
                 conf.int = TRUE,
                 cor.coef = TRUE, cor.method = "pearson",
                 xlab="Mean Coverage",
                 ylab="Portion of Genome Covered")

breadth <- ggplot(bwa,aes(x=Breadth_Coverage))+
  geom_boxplot(fill="#87BBA2", width=0.01,position= position_nudge(y=.1))+
  geom_histogram(aes(y=stat(density)), size=1, fill="#87BBA2", alpha=0.8)+
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background=element_rect(colour='black', fill='white')
  )+
  xlab("Portion of Reference Genome Covered (%)")+
  ylab("Density")

# Checking mapping qualities
mapqmelt <- gather(bwa,# having raw reads is pointless 
                   "Quality",
                   "Reads",
                   -c(Accession,Mean_Coverage,Breadth_Coverage))
mapqmelt$Quality <- factor(mapqmelt$Quality, levels=c("Mapped","MAPQ>0","MAPQ>9","MAPQ>19","MAPQ>29"))

ggplot(mapqmelt, aes(x=Quality, y=Reads,group=Quality))+
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background=element_rect(colour='black', fill='white')
  )+
  geom_boxplot()+
  xlab("Quality Bracket")+
  ylab("Read Count")

ggplot(subset(mapqmelt, Quality != "Mapped"), aes(x=Quality, y=Reads,group=Quality))+
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background=element_rect(colour='black', fill='white')
  )+
  geom_boxplot()+
  geom_point(position=position_dodge2(0.2))+
  #geom_line(aes(group=Accession)) +
  xlab("Quality Bracket")+
  ylab("Read Count")


#--------------------------------------------------------------------#
#### Quantifying Loss of Mapped Reads Following PCR Deduplication ####
#--------------------------------------------------------------------#

# Need to double check why bwa$Accession is longer than dup$Accession

dup <- read.csv(here("./Mapped Read Statistics/deduplicated_mapping_stats.csv"), header=FALSE)
names(dup) <- c("Accession","Percent_Mapped","Percent_Mapped_Paired","Alignments","MAPQ>0","MAPQ>9","MAPQ>19","MAPQ>29","Mean_Coverage","Mean_Relative_Coverage","Breadth_Coverage")

# Get the deduplicated stats in the correct order with respect to the stats from
#   the unprocessed dataset
# By running left_join, the deduplicated dataset can be reorder with respect to 
#   order of accessions in the preprocessed bwa dataset
dup_reordered <- left_join(data.frame(Accession=bwa$Accession), dup,
                    by="Accession")

# Get summary statistics for a barplot for the deduplicated data
dup_melt <- gather(subset(dup_reordered, select=c("Accession","MAPQ>0","MAPQ>9","MAPQ>19","MAPQ>29")),
                   "Quality_Bracket",
                   "Reads",
                   -Accession)

dup_summary <- dup_melt %>%
  group_by(Quality_Bracket) %>%
  summarise(Mean=mean(Reads),Sdev=sd(Reads))

# Make sure the quality brackets are set to a correctly ordered factor
dup_summary$Quality_Bracket <- factor(dup_summary$Quality_Bracket, levels=c("MAPQ>0","MAPQ>9","MAPQ>19","MAPQ>29"))

# Get a similar summary dataframe for the untreated data
bwa_melt <- gather(subset(bwa, select=c("Accession","MAPQ>0","MAPQ>9","MAPQ>19","MAPQ>29")),
                   "Quality_Bracket",
                   "Reads",
                   -Accession)

bwa_summary <- bwa_melt %>%
  group_by(Quality_Bracket) %>%
  summarise(Mean=mean(Reads),Sdev=sd(Reads))

bwa_summary$Quality_Bracket <- factor(bwa_summary$Quality_Bracket, levels=c("MAPQ>0","MAPQ>9","MAPQ>19","MAPQ>29"))

# Draw the barplot
ggplot()+
  geom_bar(data=bwa_summary,
           aes(x=Quality_Bracket, y=Mean), stat="identity", fill="#38799E")+
  geom_errorbar(data=bwa_summary,
                aes(x=Quality_Bracket, y=Mean, ymin=Mean-Sdev, ymax=Mean+Sdev))+
  geom_bar(data=dup_summary,
           aes(x=Quality_Bracket, y=Mean), stat="identity", fill="#87BBA2")+
  geom_errorbar(data=dup_summary,
                aes(x=Quality_Bracket, y=Mean, ymin=Mean-Sdev, ymax=Mean+Sdev))


### SNP statistics

snps <- read.csv(here("./SNP Statistics/SNP_log.csv"), header=F)
names(snps) <- c("Scaffold", "Position", "Chunk", "SNPs")

# get the total number of SNPs
sum(snps$SNPs)

# See the density of snps across different scaffolds
plot(snps$SNPs ~ snps$Chunk, col=factor(snps$Scaffold))

# Get a way to test for a difference in density between scaffolds
anova <- aov(SNPs ~ Scaffold, data=snps)
summary(anova)

# Look only at significant comaprisosn
pairwise <- TukeyHSD(anova, conf.level = 0.95)
subset(pairwise$Scaffold, pairwise$Scaffold[,4] < 0.05)

ggplot(data=snps)+
  geom_boxplot(aes(x=Scaffold,y=SNPs))
