# Ray-Simulation




acrylic_diameter=0.1;
acrylic_length=0.443;
dist_horn_to_acrylic=0.1 #Note: Estimate
acrylic_radius=acrylic_diameter/2;


max_angle=atan(acrylic_radius/dist_horn_to_acrylic);



n_acrylic=1.5;
n_air=1;
#max_angle=pi/2;


theta_i1=seq(from=0,to=(max_angle),by=0.0001);

incident_r=tan(theta_i1)*dist_horn_to_acrylic;
exit_r=rep(NaN,length=length(incident_r))
exit_angle=rep(NaN,length=length(incident_r))


theta_t1=n_air*sin(theta_i1)/n_acrylic;
numel_reflect=rep(0,length(theta_i1))

theta_i2=pi/2-theta_t1;

dist_horz_acrylic=matrix(0,nrow=length(theta_t1),ncol=4)
dist_total_acrylic=matrix(0,nrow=length(theta_t1),ncol=4)

dist_horz_acrylic[,1]<-(acrylic_radius-incident_r)/tan(theta_t1);
for (i in 1:length(dist_horz_acrylic[,1])){
  if (dist_horz_acrylic[i,1]>acrylic_length){
    dist_horz_acrylic[i,1]=acrylic_length
    exit_r[i]=tan(theta_t1[i])*acrylic_length+incident_r[i]
    numel_reflect[i]=numel_reflect[i]+0;
    dist_total_acrylic[i,1]=sqrt((exit_r[i]-incident_r[i])^2+(acrylic_length)^2)
  exit_angle[i]=theta_t1[i]}
    else{
    numel_reflect[i]=numel_reflect[i]+1;
    dist_total_acrylic[i,1]=sqrt((acrylic_radius-incident_r[i])^2+(dist_horz_acrylic[i,1])^2)}
  }





theta_t2=asin((n_acrylic/n_air)*sin(theta_i2));

theta_r2=theta_i2;
theta_i3=theta_r2;

dist_horz_acrylic[,2]=(2*acrylic_radius)*tan(theta_r2);

for (i in 1:length(dist_horz_acrylic[,2])){
  if(numel_reflect[i]<1){}
  else{
  if ((dist_horz_acrylic[i,1]+dist_horz_acrylic[i,2])>acrylic_length){
    
    dist_horz_acrylic[i,2]=acrylic_length-dist_horz_acrylic[i,1]
    exit_r[i]=acrylic_radius-(acrylic_length-dist_horz_acrylic[i,1])/tan(theta_r2[i]);
    
    numel_reflect[i]=numel_reflect[i]+0;
    dist_total_acrylic[i,2]=sqrt((exit_r[i]-acrylic_radius)^2+(acrylic_length-dist_horz_acrylic[i,1])^2)
  exit_angle[i]=-(pi/2-theta_r2[i])}
  else{
    numel_reflect[i]=numel_reflect[i]+1;
    dist_total_acrylic[i,2]=sqrt((2*acrylic_radius)^2+(dist_horz_acrylic[i,2])^2)}
}}

theta_r3=theta_i3;
theta_i4=theta_r3;
dist_horz_acrylic[,3]=(2*acrylic_radius)*tan(theta_r3);


for (i in 1:length(dist_horz_acrylic[,3])){
  if(numel_reflect[i]<2){}
  else{
    if ((dist_horz_acrylic[i,1]+dist_horz_acrylic[i,2]+dist_horz_acrylic[i,3])>acrylic_length){
      dist_horz_acrylic[i,3]=acrylic_length-dist_horz_acrylic[i,1]-dist_horz_acrylic[i,2]
      exit_r[i]=-acrylic_radius+(acrylic_length-dist_horz_acrylic[i,1]-dist_horz_acrylic[i,2])/tan(theta_r3[i]);
      
      numel_reflect[i]=numel_reflect[i]+0;
      dist_total_acrylic[i,3]=sqrt((exit_r[i]+acrylic_radius)^2+(acrylic_length-dist_horz_acrylic[i,1]-dist_horz_acrylic[i,2])^2)
    exit_angle[i]=pi/2-theta_r3[i]}
    else{
      numel_reflect[i]=numel_reflect[i]+1;
      dist_total_acrylic[i,2]=sqrt((2*acrylic_radius)^2+(dist_horz_acrylic[i,2])^2)}
  }}


dist_horz_acrylic_tot=dist_horz_acrylic[,1]+dist_horz_acrylic[,2]+dist_horz_acrylic[,3]
dist_tot=dist_total_acrylic[,1]+dist_total_acrylic[,2]+dist_total_acrylic[,3]
out_angle=asin(sin(exit_angle)*n_acrylic/n_air)

dist_out=0.15 #Note: Estimate

pos_out=tan(out_angle)*dist_out

hist_pos=hist(pos_out,15)
plot(theta_i1,dist_tot,xlab='Incident Angle (Radians)', ylab='Total Distance Traveled within Acrylic Waveguide (m)',main='Distance Travelled Within Acrylic Waveguide vs. Incident Angle')

plot(theta_i1,numel_reflect,xlab='Incident Angle (Radians)', ylab='Number of Internal Reflections',main='Number of Internal Reflections vs. Incident Angle')
