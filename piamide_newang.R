library('bio3d','dplyr', 'ggplot2')
library('pracma') # for the the cross function of normalizing the vector
rm(list=ls())
args <- commandArgs(trailingOnly=TRUE)
setwd(args[1])
temp <- args[2]; print(paste0('this is the temp environment: ', temp))
test <- args[3]; print(paste0('this is the test number: ', test))
print('need to load in a topology and trajectory that has already had the hydrogens removed')
message <- (paste('easily done in cpptraj:','strip @H=',
        'trajout xxx.pdb onlyframes 1',
        'trajout xxxx.dcd', sep = '\n'))
cat(message)
tpr <- read.pdb('tpr-noh.pdb')
trj <- read.dcd('trj-noh.dcd')
#noh <- trim.pdb(tpr, atom.select(tpr, string='noh'))
#sel <- seq(0,10000, by = 10)
#trj_small <- trim(trj,row.inds=sel)

mut.inds <- atom.select(tpr, string='noh')

xyz = fit.xyz(fixed = tpr$xyz, mobile = trj,
              fixed.inds = mut.inds$xyz, mobile.inds = mut.inds$xyz)

info <- 'mass of carbon in ring = 12.0107'
cat(info)
info2 <- 'masses in CON of Y17/F18: 12.0107, 15.9994, 14.0067, respectively'
cat(info2)
oxy <- 15.9994; nitro <- 14.0067; carb <- 12.0107

distance_cutoff <- 4.5
angle_cutoff <- 140

data <- data.frame(frames = 1:length(xyz[,1]), double = FALSE, in_dist = 0, dist = 0, 
                   in_ang =0, angle = 0, in_both_ang = 0,
                   ring_pos = 0,amide_pos = 0, normed_pos = 0)

#ring.inds <- c(75,76,77,78,79,80) # extract the ring indices for N13F
ring.inds <- c(79,80,81,82,83,84)
ringname <- c('CG', 'CD1', 'CE1','CZ','CE2','CD2')
ring <- atom.select(tpr, elety = ringname, resno = 12) #limit the selection these specific atoms
ring.plane <- atom.select(tpr,eleno=c(79, 84, 83)) # obtain the plane

#amide.inds <- c(111,112,113) # extract the atom indice of the amide region Y17/F18
amide.inds <- c(118,119,120)
amide <- atom.select(tpr, eleno = amide.inds)
#ringxyz <- tpr$xyz[,ring$xyz] # tpr$xyz[,ring$xyz] if using read.pdb file, tpr
#ring_pos <- com.xyz(ringxyz, mass = c(carb, carb, carb, carb, carb, carb))


message2 <- paste('need to get the normal plane of the amide',
                  'we want a position above the com of the amide',
                  'then we can go down to the amide com then over to ring com',
                  'with the angle at the amide com', sep = '\n')


#amidexyz <- tpr$xyz[,amide$xyz]  # tpr$xyz[,amide$xyz] if using read.pdb file, tpr 
#amide_pos <- com.xyz(amidexyz, mass = c(carb, oxy, nitro))


for ( frame in 1:length(xyz[,1])) {
  top <- xyz[frame,]
  ringxyz <- top[ring$xyz] # tpr$xyz[,ring$xyz] if using read.pdb file, tpr
  ring_pos <- com.xyz(ringxyz, mass = c(carb, carb, carb, carb, carb, carb))
  
  cg <- ringxyz[1:3] # extract the coordinates of the CG atom in the ring for angle measurements
  
  
  amidexyz <- top[amide$xyz] # tpr$xyz[,amide$xyz] if using read.pdb file, tpr 
  amide_pos <- com.xyz(amidexyz, mass = c(carb, oxy, nitro))
  data$amide_pos[frame] <- paste(amide_pos[1],amide_pos[2],amide_pos[3],sep = ' ')
  # assemble the distance matrix between the ring and amide: 
  coords_matrix <- rbind(ring_pos, amide_pos)
  distance <- round(as.numeric(dist(coords_matrix)),4)
  # save these to the dataframe at the current frame -- also save position
  data$in_dist[frame] <- ifelse(distance <= distance_cutoff, TRUE, FALSE)
  data$dist[frame] <- distance
  data$ring_pos[frame] <- paste(ring_pos[1],ring_pos[2],ring_pos[3],sep = ' ')
  # calculate the angle by getting the normal plane
  amc <- com.xyz(amidexyz[1:3], mass = carb); amc # extract the coordinates of the C atom in the amide for angle measurements
  amo <- com.xyz(amidexyz[4:6], mass = oxy); amo
  amn <- com.xyz(amidexyz[7:9], mass = nitro); amn
  
  # vector subtraction to get the difference in 3D space  
  amo_amc <- amo - amc; amo_amc
  amo_amn <- amo - amn ; amo_amn
  #get the crossproduct and normalize the perpendicular vector 
  normal_vector <- round(cross(amo_amc, amo_amn), 3); normal_vector
  normed <- round(normal_vector / sqrt(sum(normal_vector^2)),4); normed
  #finally, displace amide COM by the normalized vector
  ref <- round(amide_pos + normed, 3); ref
  #print('getting ring plane')
  ringplane.xyz <- top[ring.plane$xyz]
  cg <- com.xyz(ringplane.xyz[1:3],mass=carb); cg
  cz <- com.xyz(ringplane.xyz[4:6],mass=carb); cz
  ce2 <- com.xyz(ringplane.xyz[7:9],mass=carb); ce2
  cz_cg <- cz - cg
  cz_ce2 <- cz - ce2
  normal_vector.ring <- round(pracma::cross(cz_cg,cz_ce2), 3); #normal_vector
  normed.ring <- round(normal_vector.ring / sqrt(sum(normal_vector.ring^2)),4); #normed
  #finally, displace amide COM by the normalized vector
  ref.ring <- round(ring_pos - normed.ring, 3); #ref
  angle.ringplane <- angle.xyz(c(ring_pos,ref.ring,amide_pos))
  #if (angle.ringplane < 90) {
   # ref.ring <- round(ring_pos + normed.ring,3)
    #angle.ringplane <- angle.xyz(c(ring_pos,ref.ring,amide_pos))
  #}
  data$normed_pos[frame] <- paste(ref.ring[1],ref.ring[2],ref.ring[3],sep = ' ') #was list(ref)
  #calculate angle: normed_vec, amide_com, ring_com
  angle <- angle.xyz(c(ref, amide_pos, ring_pos)); angle
  data$angle[frame] <- angle.ringplane
  data$in_ang[frame] <- ifelse(angle.ringplane >= angle_cutoff | angle.ringplane <= 40,
                               TRUE, FALSE)
  data$double[frame] <- ifelse(data$in_ang[frame] == TRUE & distance <= distance_cutoff, 
                               TRUE, FALSE)
} 
angle_dist_true <- length(subset(data, double == TRUE)[,1])/length(xyz[,1])*100
print(paste0(angle_dist_true, '% of the data is represented within the angle and distance cutoffs of: ', 
             angle_cutoff,' degrees and ', distance_cutoff, ' angstroms'))
confirm <- subset(data, double ==TRUE)

write.table(confirm,paste0(temp,'-',test,'-','y17f18-piamide-data_CONFIRM.xvg'), sep = '\t')
write.table(data,paste0(temp,'-',test,'-','y17f18-piamide-data.xvg'), sep = '\t')

print('making table for the pi-amide interaction')
angle.avg <- round(mean(data$angle), 3); angle.avg
angle.sd <- round(sd(data$angle),3); angle.sd
dist.avg <- round(mean(data$dist), 3); dist.avg
dist.sd <- round(sd(data$dist),3); dist.sd
perc.cut <- paste0(angle_dist_true,'%'); perc.cut
data_out <- data.frame('Test' = paste0(temp,'-',test),'Percent in Cutoff' =perc.cut, 'Average Dist'=dist.avg, 'Average Angle' = angle.avg,
                       'Std. Deviation of Distance' = dist.sd, 'Std. Deviation of Angle' = angle.sd)
print('here is the final summary of this data: ')
data_out
write.table(data_out, paste0('summary_pi-amide_',temp,'-',test,'.xvg'), sep = '\t')

fig_title <- 'N13F ring and Y17/F18 C,O, and N meeting groups'

library('ggplot2')
dist_fig <- ggplot2::ggplot(data, aes(x=frames, color=protein)) +
 # geom_line(aes(y=X11P, color="X11P"), lwd=1) +
  geom_line(aes(y=dist, color = 'pi-amide'), lwd=1) +
  labs(x="Frame (1000 = 10ns)", y = "Distance (Å)", title=paste('Distance of ', fig_title, sep = '')) +
  theme_classic() + theme(legend.position = 'top', legend.title = element_blank(),
                          legend.background = element_rect(fill = NA, color= NA),
                          legend.key = element_rect(fill = "transparent", color= NA),
                          legend.text = element_text(size = 18), legend.key.size = unit(1,'cm'), 
                          axis.text = element_text(face='bold', size=20, color = 'black'),
                          axis.title = element_text(face='bold', size = 20),
                          plot.title = element_text(face='bold', hjust=0.5, size = 20),
                          panel.border = element_rect(fill = NA, linewidth=1)) +
  scale_x_continuous(breaks = seq(0,length(data[,1]-1), by = length(data[,1])/10)) +
  scale_y_continuous(breaks=seq(0,10,by=2), limits = c(0,10),
                     sec.axis = dup_axis(name = NULL))
dist_fig

check <- which(data$angle < 10 & data$dist <= 4.0 ); check
data$angle[2]
write.pdb(pdb = tpr, xyz = xyz[999,], file = 'test.pdb')

ggsave(paste0(temp,'-',test,'-','pnq-y17-f18-piamide-distance.png'), dist_fig, width=10, height = 8)


ang_fig <- ggplot2::ggplot(data, aes(x=frames, color=protein)) +
  # geom_line(aes(y=X11P, color="X11P"), lwd=1) +
  geom_line(aes(y=angle, color = 'pi-amide'), lwd=1) +
  labs(x="Frame (1000 = 10ns)", y = "Angle (Degrees)", title=paste('Angle of ', fig_title, sep = '')) +
  theme_classic() + theme(legend.position = 'top', legend.title = element_blank(),
                          legend.background = element_rect(fill = NA, color= NA),
                          legend.key = element_rect(fill = "transparent", color= NA),
                          legend.text = element_text(size = 18), legend.key.size = unit(1,'cm'), 
                          axis.text = element_text(face='bold', size=20, color = 'black'),
                          axis.title = element_text(face='bold', size = 20),
                          plot.title = element_text(face='bold', hjust=0.5, size = 20),
                          panel.border = element_rect(fill = NA, linewidth=1)) +
  scale_x_continuous(breaks = seq(0,length(data[,1]-1), by = length(data[,1])/10)) +
  scale_y_continuous(breaks=seq(0,180,by=15), limits = c(0,180),
                     sec.axis = dup_axis(name = NULL))
ang_fig


ggsave(paste0(temp,'-',test,'-','pnq-y17-f18-piamide-angle-343_1.png'), ang_fig, width=10, height = 8)

