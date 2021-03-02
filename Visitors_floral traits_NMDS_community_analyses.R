
#R script used to run NMDS analyses of insect visitor community fit to floral trait vectors, in order to identify most important floral traits 
#for modeling visitation rates (see Multivariate_visitation_rate_models.R) 
#In: Drager et al. 'What structures diurnal visitation rates to flowering trees in an Afrotropical lowland rain-forest understory?' 

library(vegan)
library(ggplot2)
library(dplyr)

#Visitors

read.csv("KUT_visitor_comm_prop.csv", header=T, row.names=1,stringsAsFactors = FALSE)->visp
as.matrix(visp)->visMp
t(visMp)->tvisMp

p2_NMDS=metaMDS(tvisMp,k=2, try=1000, trymax=1000, autotransform=T, wascores=T)  #Stress: 0.1029327 #2D, when using prop matrix 
p3_NMDS=metaMDS(tvisMp,k=3, try=1000, trymax=1000, autotransform=T, wascores=T)  #Stress: 0.05161172 #3D, when using prop matrix

stressplot(p3_NMDS)  #all points along the line --> 3D works well, use this one

site.scores.3 <- as.data.frame(scores(p3_NMDS)) 
site.scores.3$site <- rownames(site.scores.3)
species.scores.3 <- as.data.frame(scores(p3_NMDS, "species")) 
species.scores.3$species <- rownames(species.scores.3)

#Floral Traits

read.csv("KUT_floral_traits.csv", header=T, row.names=1,stringsAsFactors = FALSE)->traits
as.matrix(traits)->tM 
t(tM)->ttM  #transpose

#fitting traits onto visitor group-tree spp ordination.

vf12 <- envfit(p3_NMDS,ttM,choices=1:2, permu = 9999)  #two significant traits (sweet,fermented) --> choose this for figure
vf13 <- envfit(p3_NMDS,ttM,choices=1:3, permu = 9999)  #one marginaly significant traits (fermented)
vf23 <- envfit(p3_NMDS,ttM,choices=2:3, permu = 9999)  #no significant traits


#moving 2D plot to ggplot:

#organize trait data:
trait.NMDS.data<-as.data.frame(ttM)

trait.NMDS.data$NMDS1<-p3_NMDS$points[ ,1] #this puts NMDS scores for the plots into a new dataframe. 
trait.NMDS.data$NMDS2<-p3_NMDS$points[ ,2] 

vis.NMDS.data<-as.data.frame(tvisMp)
stems<-colSums(vis.NMDS.data)
spps <- data.frame(scores(p3_NMDS, display = "species")) #dataframe of species scoes for plotting
spps$species <- row.names(spps) # making a column with species names
spps <- cbind(spps, data.frame(GROUP = c("ANTS","BEES","BEETLES","BUGS","DELICATE_FLIES","ROBUST_FLIES")))
spps$colsums <- stems #adding the colSums from above
spps<-spps[!is.na(spps$NMDS1) & !is.na(spps$NMDS2),] #removes NAs
spps$GROUP <-as.character(spps$GROUP)

# data for the envfit arrows:

trait.scores <- as.data.frame(scores(vf12n, display = "vectors")) #extracts relevant scores from envifit
#trait.scores<- cbind(trait.scores, trait.variables = rownames(trait.scores)) #and then gives them their names
trait.scores<-mutate(trait.scores, floral_trait=c("morphology", "scent","scent","scent","scent","scent","scent","color", "color","color","color","color","color", "color", "visible nectar"))
trait.scores<-mutate(trait.scores, trait=c("restrictive", "sweet","spicy","fermented","carrion","fruity","unscented","yellow", "red","white","orange","green","pink", "maroon", "visible nectar"))

# finally plotting: 

Trait_palette1 = c("morphology"="#56B4E9", "scent"= "#009E73", "scent"="#009E73", "scent"="#009E73", "scent"="#009E73", "scent"="#009E73", "color"="#F0E442", "color"= "#F0E442","color"= "#F0E442","color"="#F0E442", "color"= "#F0E442","color"="#F0E442", "color"="#F0E442", "color"="#009E73","visible nectar"="#CC79A7")
scales::show_col(Trait_palette1)

mult <- 2 #multiplier for the arrows and text for envfit below. You can change this and then rerun the plot command.
(korup.nmds.gg1 <- ggplot(data = trait.NMDS.data, aes(y = NMDS2, x = NMDS1))+ #sets up the plot. brackets around the entire thing to make it draw automatically
    geom_point(size = 3,color="darkgray") + #puts the site points in from the ordination, shape determined by site, size refers to size of point
    geom_text(data=spps, aes(x=spps$NMDS1, y=spps$NMDS2, label=GROUP), size = 4,color="black")+ #labelling the species. hjust used to shift them slightly from their points for readability.
    geom_segment(data = trait.scores,
                 aes(x = 0, xend = mult*NMDS1, y = 0, yend = mult*NMDS2, color=floral_trait),
                 arrow = arrow(length = unit(0.25, "cm"))) + 
    scale_color_manual(values=Trait_palette1)+ #arrows for envfit.  doubled the length for similarity to the plot() function. NB check ?envfit regarding arrow length if not familiar with lengths
    geom_text(data = trait.scores, #labels the environmental variable arrows * "mult" as for the arrows
              aes(x = mult*NMDS1, y = mult*NMDS2, label=trait),
              size = 4,
              position=position_jitter(width=.04, height=.02)))+
  annotate("text",x =-1.12 ,y = 0.5,label="*", size=8) + 
  annotate("text",x =1.20 ,y = -0.08,label="*", size=8) + 
  coord_cartesian(xlim = c(-1.2,1.2), ylim= c(-0.75,1.1))+  # NB this changes the visible area of the plot only (this is a good thing, apparently). Can also specify ylim. Here in case you want to set xaxis manually.
  theme_bw()



t(tMn)->ttMn