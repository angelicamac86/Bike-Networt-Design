# Bike-Networt-Design for Flat cities
Models to bike network design
% TRAVEL DEMAND DISTRIBUTIONS
x=-20:1:19
y=-20:1:19
Distrib0=zeros (40,40)
Distrib1(1:40,1:40) =1.0
for i=1:40
    for j=1:40
       Distrib2(i, j) =max (1,4-x (i). *x (i)/80-y (j). *y (j)/80); 
    end
end

 Ao=Distrib1; % origins users A
 Ad=Distrib1; % destinations users A
 Bo=Distrib1; % origins users B
 Bd=Distrib1; % destinations users B
 Co=Distrib1; % origins users C
 Cd=Distrib1; % destinations users C

%INPUTS
     Dx = 15;
     Dy = 8;
     alpha = 0.7;
     deltapt = 7.11; 
     vw = 4;
     vb1 = 12;
     vb2 = 15;
     Cl = 1.24; %Cl=1.24
     Cs = 0.46; %Cs = 0.46;
     Cm = 0.041; %Cm = 0.041;
     Crepoteams=5.561; %Coste de cada equipo de repositioning. En BCN es 22.8 Euros/h-team Ratio 11.40/ 2.78
     Coper=0.15534; % Coste de operación (mantenimiento) de cada viaje en bici compartida
     vancap=32; % Capacidad de bicicletas en la furgoneta de repositioning 
     factordis=0.339;
     vvan=20.6; %velocidad de la furgoneta
     deltaload=37.5/3600 %tiempo de carga/descarga de una bicicleta
     Vot = 2.78;
     taux=0/3600; 
     horizon=3;
     tauy=0/3600;
     tetha=0.25; 
     deltax = Dx / 40;
     Deltay = Dy / 40;
     C=2100; 
     ixe=1; %si ixe=1 hay backtracking     
     
% SERVICE AREA DEFINITION

     dx_in = Dx * alpha;
     dy_in = Dy * alpha;
     nx_in = floor (dx_in / deltax);
     ny_in = floor (dy_in /Deltay);
     nx_out = floor ((40 - nx_in) / 2);
     ny_out = floor ((40 - ny_in) / 2);
     kx=0;
     if nx_in == ceil (nx_in / 2) *2
        kx = 0;
     else
        kx = 1;
     end 
     if ny_in == ceil (ny_in / 2) * 2 
        ky = 0;
     else
       ky = 1;
     end 
     ax = nx_out * deltax + kx * 0.5 * deltax;
     bx = (nx_in + nx_out) * deltax + kx * 0.5 * deltax;
     ay = ny_out * Deltay + ky * 0.5 * Deltay;
     by = (ny_in + ny_out) * Deltay + ky * 0.5 * Deltay; 
     
% INITIALIZATION OF TERMS THAT DO NOT DEPEND ON THE DECISION VARIABLES
      g_aVxpos=zeros (40,40);
      g_aVxneg=zeros (40,40);
      g_bVxpos=zeros (40,40);
      g_bVxneg=zeros (40,40);
      g_aHxpos=zeros (40,40);
      g_aHxneg=zeros (40,40);
      g_bHxpos=zeros (40,40);
      g_bHxneg=zeros (40,40);
      g_aHypos=zeros (40,40);
      g_aHyneg=zeros (40,40);
      g_bHypos=zeros (40,40);
      g_bHyneg=zeros (40,40);
      g_aVypos=zeros (40,40);
      g_aVyneg=zeros (40,40);
      g_bVypos=zeros (40,40);
      g_bVyneg=zeros (40,40);
      q_axpos=zeros (40,40);
      q_aypos=zeros (40,40);
      q_axneg=zeros (40,40);
      q_ayneg=zeros (40,40);
      q_bxpos=zeros (40,40);
      q_bypos=zeros (40,40);
      q_bxneg=zeros (40,40);
      q_byneg=zeros (40,40);
      fc12=zeros (40,40);
      fc_primer=zeros (40,40);
      fc_tercer=zeros (40,40);
      fc3=zeros (40,40);
      fc_segundo=zeros (40,40);
      fc_cuarto=zeros (40,40);
      fa= Ad. *sum (Ao, 'all') +Ao. *sum (Ad,'all'); 
      fa=fa. *deltax*Deltay;
      fa;
      fb=Bd. *sum (Bo, 'all') +Bo. *sum (Bd, 'all');
      fb=fb. *deltax*Deltay;
      fb;
      fb_out= Bo. *sum (Bd, 'all');
      fb_out=fb_out. *deltax*Deltay
      fb_in= Bd. *sum (Bo, 'all')
      fb_in=fb_in. *deltax*Deltay
      
% FLOW CALCULATIONS
  
  for i=1:40
    for j=1:40    
  
              for k1 = 1: j
               for k2 = j: 40
                 if k1~=k2
                   for k3 = 1: i
                    for k4 = i: 40
                      if k3~= k4 
                      g_aVxpos (i, j) = g_aVxpos (i, j) + Ao (k3,k1) *Ad(k4, k2) * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / Deltay + Ao (k4,k1) *Ad(k3,k2) * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / Deltay;
                       
                        g_aVxneg (i, j) = g_aVxneg (i, j) + Ao (k3,k2) *Ad(k4, k1) * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / Deltay + Ao (k4,k2) *Ad(k3,k1)  * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / Deltay;
                     
                        g_bVxpos (i, j) = g_bVxpos (i, j) + Bo (k3,k1) *Bd(k4, k2) * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / Deltay + Bo (k4,k1) *Bd(k3,k2) * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / Deltay;
                       
                        g_bVxneg (i, j) = g_bVxneg (i, j) + Bo (k3,k2) *Bd(k4, k1) * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / Deltay + Bo (k4,k2) *Bd(k3,k1) * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / Deltay;
                      end
                    end
                  end
                  
     for k5 = 1: 40
                    g_aHxpos (i, j) = g_aHxpos (i, j) + Ao (i,k1) *Ad(k5, k2) * deltax * deltax * Deltay / 2 / (k2 - k1) * (k2 - j) + Ao (k5,k1) *Ad(i, k2) * deltax * deltax * Deltay / 2 / (k2 - k1) * (j - k1);
                       
                      g_aHxneg (i, j) = g_aHxneg (i, j) + Ao (i,k2) *Ad(k5, k1) * deltax * deltax * Deltay / 2 / (k2 - k1) * (j - k1) + Ao (k5,k2) *Ad(i, k1) * deltax * deltax * Deltay / 2 / (k2 - k1) * (k2 - j);

                      g_bHxpos (i, j) = g_bHxpos (i, j) + Bo (i,k1) *Bd(k5, k2) * deltax * deltax * Deltay / 2 / (k2 - k1) * (k2 - j) + Bo (k5,k1) *Bd(i, k2)  * deltax * deltax * Deltay / 2 / (k2 - k1) * (j - k1);
                       
                   g_bHxneg (i, j) = g_bHxneg (i, j) + Bo (i,k2) *Bd(k5, k1) * deltax * deltax * Deltay / 2 / (k2 - k1) * (j - k1) + Bo (k5,k2) *Bd(i,k1) * deltax * deltax * Deltay / 2 / (k2 - k1) * (k2 - j);
                  end
                 end
               end
             end
            
  for k1 = 1: i
               for k2 = i: 40
                 if k1~=k2
                   for k3 = 1: j
                    for k4 = j: 40
                      if k3~=k4
                     g_aHypos (i, j) = g_aHypos (i, j) + Ao (k1,k3) *Ad(k2, k4)  * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / deltax + Ao (k1,k4) *Ad(k2, k3) * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / deltax;
                       
                       g_aHyneg (i, j) = g_aHyneg (i, j) + Ao (k2,k4) *Ad(k1, k3)  * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / deltax + Ao (k1,k3) *Ad(k2, k4) * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / deltax;
                      
                       g_bHypos (i, j) = g_bHypos (i, j) + Bo (k1,k3) *Bd(k2, k4) * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / deltax + Bo (k1,k4) *Bd(k2, k3) * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / deltax;
                      
                       g_bHyneg (i, j) = g_bHyneg (i, j) + Bo (k2,k4) *Bd(k1, k3) * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / deltax + Bo (k1,k3) *Bd(k2, k4) * deltax * deltax * Deltay * Deltay / 2 / (k4 - k3) / deltax;
                      end
                    end
                  end

                  for k5 = 1: 40
                      g_aVypos (i, j) = g_aVypos (i, j) + Ao (k1,j) *Ad(k2, k5)  * deltax * Deltay * Deltay / 2 / (k2 - k1) * (k2 - i) + Ao (k1,k5) *Ad( k2,j) * deltax * Deltay * Deltay / 2 / (k2 - k1) * (i - k1);

                      g_aVyneg (i, j) = g_aVyneg (i, j) + Ao (k2,j) *Ad(k1, k5) * deltax * Deltay * Deltay / 2 / (k2 - k1) * (i - k1) + Ao (k2,k5) *Ad( k1,j) * deltax * Deltay * Deltay / 2 / (k2 - k1) * (k2 - i);

                      g_bVypos (i, j) = g_bVypos (i, j) + Bo (k1,j) *Bd(k2, k5) * deltax * Deltay * Deltay / 2 / (k2 - k1) * (k2 - i) + Bo (k1,k5) *Bd( k2,j) * deltax * Deltay * Deltay / 2 / (k2 - k1) * (i - k1);

                      g_bVyneg (i, j) = g_bVyneg (i, j) + Bo (k2,j) *Bd(k1, k5) * deltax * Deltay * Deltay / 2 / (k2 - k1) * (i - k1) + Bo (k2,k5) *Bd(k1,j) * deltax * Deltay * Deltay / 2 / (k2 - k1) * (k2 - i);
                  end
                 end
               end
             end

         % user C
           if (i<=ny_out)|(i>=ny_out + 1 + ny_in)|(j<=nx_out)|(j>=nx_out + 1 + nx_in)
             for k1 = nx_out: (nx_in + nx_out)
              for k2 = ny_out: (ny_in + ny_out)
                 fc_primer (i, j) = fc_primer (i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay;
                 fc_tercer (i, j) = fc_tercer (i, j) + Co (k2,k1) *Cd(i, j) * deltax * Deltay;
                 fc12(i, j) = fc12(i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay + Co (k2,k1) *Cd(i, j) * deltax * Deltay;
                end
              end
            end
             
             if (i<=ny_out) 
                 if (j<=nx_out) 
                    zone=1;
                 elseif (j>=nx_out + 1 + nx_in)
                    zone=3;
                 else
                    zone=2;  
                 end   
             end

             if (i>=ny_out + 1 + ny_in) 
                 if (j<=nx_out) 
                    zone=7;
                 elseif (j>=nx_out + 1 + nx_in)
                    zone=9;
                 else
                    zone=8;
                 end   
             end

             if ((i>ny_out) & (i<=ny_in + ny_out)) 
                 if (j<=nx_out) 
                    zone=4;
                 elseif (j>=nx_out + 1 + nx_in)
                    zone=6;
                 else
                    zone=5;
                 end   
             end                    
        %Dibujo 2.1
             if (zone==1)|(zone==2)|(zone==3)|(zone==4)|(zone==6)|(zone==7)|(zone==8)|(zone==9) 
              %destino zona=1   
              for k1 = 1: nx_out
               for k2 = 1: ny_out
                   fc3(i, j) = fc3(i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay + Co (k2,k1) *Cd(i, j) * deltax * Deltay;
                   fc_segundo (i, j) = fc_segundo (i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay;
                   fc_cuarto (i, j) = fc_cuarto (i, j) + Co (k2,k1) *Cd(i, j)  * deltax * Deltay;
               end
              end
             end

            if (zone==1)|(zone==2)|(zone==3)|(zone==4)|(zone==6)|(zone==7)|(zone==8)|(zone==9)
               %destino zona=7   
              for k1 = 1: nx_out
               for k2 = ny_out + 1 + ny_in: 40
                   fc3(i, j) = fc3(i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay + Co (k2,k1) *Cd(i, j) * deltax * Deltay;
                   fc_segundo (i, j) = fc_segundo (i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay;
                   fc_cuarto (i, j) = fc_cuarto (i, j) + Co (k2,k1) *Cd(i, j)  * deltax * Deltay;
               end
              end
             end

             if (zone==1)|(zone==2)|(zone==3)|(zone==4)|(zone==6)|(zone==7)|(zone==8)|(zone==9)
                 %destino zona=4   
              for k1 = 1: nx_out
               for k2= (ny_out + 1) : (ny_out+ ny_in) 
                   fc3(i, j) = fc3(i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay + Co (k2,k1) *Cd(i, j) * deltax * Deltay;
                   fc_segundo (i, j) = fc_segundo (i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay;
                   fc_cuarto (i, j) =fc_cuarto (i, j) + Co (k2,k1) *Cd(i, j) * deltax * Deltay;
               end
              end
             end

             if (zone==1)|(zone==2)|(zone==3)|(zone==4)|(zone==6)|(zone==7)|(zone==8)|(zone==9)
                 %destino zona=2   
              for k1 = (nx_in + 1): (nx_out + nx_in) 
               for k2 = 1: ny_out 
                   fc3(i, j) = fc3(i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay + Co (k2,k1) *Cd(i, j) * deltax * Deltay;
                   fc_segundo (i, j) = fc_segundo (i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay;
                   fc_cuarto (i, j) = fc_cuarto (i, j) + Co (k2,k1) *Cd(i, j)  * deltax * Deltay;
               end
              end
             end

             if (zone==1)|(zone==2)|(zone==3)|(zone==4)|(zone==6)|(zone==7)|(zone==8)|(zone==9)
                 %destino zona=8   
              for k1 = (nx_in + 1) : (nx_out+ nx_in) 
               for k2 = ny_out + 1 + ny_in : 40
                   fc3(i, j) = fc3(i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay + Co (k2,k1) *Cd(i, j) * deltax * Deltay;
                   fc_segundo (i, j) = fc_segundo (i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay;
                   fc_cuarto (i, j) = fc_cuarto (i, j) + Co (k2,k1) *Cd(i, j)  * deltax * Deltay;
               end
              end
             end

             if (zone==1)|(zone==2)|(zone==3)|(zone==4)|(zone==6)|(zone==7)|(zone==8)|(zone==9)
                 %destino zona=3 
              for k1 = (nx_out+ nx_in+1):40
               for k2 = 1: ny_out
                   fc3(i, j) = fc3(i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay + Co (k2,k1) *Cd(i, j) * deltax * Deltay;
                   fc_segundo (i, j) = fc_segundo (i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay;
                   fc_cuarto (i, j) = fc_cuarto (i, j) + Co (k2,k1) *Cd(i, j)  * deltax * Deltay;
               end
              end
             end
             if (zone==1)|(zone==2)|(zone==3)|(zone==4)|(zone==6)|(zone==7)|(zone==8)|(zone==9)
                 %destino zona=6 
              for k1 = (nx_out+ nx_in+1):40
               for k2 = (ny_out + 1): (ny_out+ ny_in) 
                   fc3(i, j) = fc3(i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay + Co (k2,k1) *Cd(i, j) * deltax * Deltay;
                   fc_segundo (i, j) = fc_segundo (i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay;
                   fc_cuarto (i, j) = fc_cuarto (i, j) + Co (k2,k1) *Cd(i, j)  * deltax * Deltay;
               end
              end
             end

             if (zone==1)|(zone==2)|(zone==3)|(zone==4)|(zone==6)|(zone==7)|(zone==8)|(zone==9)
                 %destino zona=9
              for k1 = (nx_out+ nx_in +1):40
               for k2 = (ny_out+ ny_in +1):40
                   fc3(i, j) = fc3(i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay + Co (k2,k1) *Cd(i, j) * deltax * Deltay;
                   fc_segundo (i, j) = fc_segundo (i, j) + Co (i, j) *Cd(k2, k1) * deltax * Deltay;
                   fc_cuarto (i, j) = fc_cuarto (i, j) + Co (k2,k1) *Cd(i, j)  * deltax * Deltay;
            end
           end
          end
        end
     end  
   svar = sqrt (1 /2/ deltapt); %svar=PT spacing/2
     
     Adem= Ao. *sum (Ad, 'all') *deltax*Deltay*deltax*Deltay;
     Bdem=Bo. *sum (Bd, 'all') *deltax*Deltay*deltax*Deltay;
     Cdem=fc3. *deltax*Deltay+fc12. *deltax*Deltay;
     Apax=sum (Adem,'all') %numero de viajes de tipo A
     Bpax=sum (Bdem,'all') %numero de viajes de tipo B
     Cpax=sum (Cdem,'all') %numero de viajes de tipo C
     
% OBJECTIVE FUNCTION CALCULATION

  for i=1:40
      i
  for j=1:40
   kkk=0;  

  % CALCULATION OF THE OBJECTIVE FUNCTION TERMS (depend on the discretization of ss, sb)
  Zmin=100000000;
   for mm2=1:40    
    xsol (2) =mm2*0.1; 
    for mm3=1:40
     xsol (3) =mm3*0.1; 
     for mm4=1:40
      xsol (4) =mm4*0.1; 
      Ns=0;
      L=0;
      A_ab=0;
      IVTT_ax=0;
      IVTT_ay=0;
      A_bw=0;
      A_bb=0;
      IVTT_bx =0;
      IVTT_by=0; 
      A_cwT=0;
      A_cb1=0;
      IVTT_c=0;
      Ns = Ns + 2 / xsol (2) / xsol (2) * deltax * Deltay;
      L = L + 2 * deltax / xsol (3) * Deltay;
      L = L + 2 * Deltay / xsol (4) * deltax;
      
      if ixe==1 % if ixe= 0 then they minimize total distance (they go direct without backtraking)
         A_ab = A_ab + fa (i, j) * (xsol (3) + xsol (4)) / 8 *(1/vb1+1/(4*vb2)) * deltax * Deltay; 
      end
      if ixe==0 % If ixe=1 then minimize walking and slow bike distance (no backtracking)
         A_ab = A_ab + fa (i, j) * (xsol (3) + xsol (4)) / 2 / vb1 * deltax * Deltay;
      end
      
      IVTT_ax = IVTT_ax + (g_aVxpos (i, j) + g_aHxpos (i, j) + g_aVxneg (i, j) + g_aHxneg (i, j)) * (1 / vb2 + taux / xsol (3)) * deltax * Deltay;
       IVTT_ay = IVTT_ay + (g_aVypos (i, j) + g_aHypos (i, j) + g_aVyneg (i, j) + g_aHyneg (i, j)) * (1 / vb2 + tauy / xsol (4)) * deltax * Deltay;
      q_axpos (i, j) = (g_aVxpos (i, j) + g_aHxpos (i, j)) * xsol (4);
      q_axneg (i, j) = (g_aVxneg (i, j) + g_aHxneg (i, j)) * xsol (4);
      q_aypos (i, j) = (g_aVypos (i, j) + g_aHypos (i, j)) * xsol (3);
      q_ayneg (i, j) = (g_aVyneg (i, j) + g_aHyneg (i, j)) * xsol (3);
      
      if ixe==0 
       A_bw = A_bw + fb (i, j) *(3*xsol (2)/4)/ vw * deltax * Deltay; 
       A_bb = A_bb + fb (i, j) * 2*(xsol (3) + xsol (4)) / 4 / vb1 * deltax * Deltay;
      end
      
      if ixe==1 
       A_bw = A_bw + fb (i, j) *(xsol (2)/3 / vw) * deltax * Deltay; 
       A_bb = A_bb + fb (i, j) * (xsol (3) + xsol (4)) * (1/ 8 / vb1+1/32/vb2) * deltax * Deltay; 
      end
      
      IVTT_bx = IVTT_bx + (g_bVxpos (i, j) + g_bHxpos (i, j) + g_bVxneg (i, j) + g_bHxneg (i, j)) * (1 / vb2 + taux / xsol (3)) * deltax * Deltay;
       IVTT_by = IVTT_by + (g_bVypos (i, j) + g_bHypos (i, j) + g_bVyneg (i, j) + g_bHyneg (i, j)) * (1 / vb2 + tauy / xsol (4)) * deltax * Deltay;
      q_bxpos (i, j) = (g_bVxpos (i, j) + g_bHxpos (i, j)) * xsol (4);
      q_bxneg (i, j) = (g_bVxneg (i, j) + g_bHxneg (i, j)) * xsol (4);
      q_bypos (i, j) = (g_bVypos (i, j) + g_bHypos (i, j)) * xsol (3);
      q_byneg (i, j) = (g_bVyneg (i, j) + g_bHyneg (i, j)) * xsol (3);
       
      % USER C
      if ixe==0
       dc1 = (3*xsol (2) / 4 *(1 / vw - 1 / vb1) + (tetha / Vot)/ (1 / vw - 1 / vb1)) ;
       dc2 = dc1+(3*xsol (2) / 4 * (1 / vb1 - 1 / vb2) + (xsol (3) +xsol (4))/4*(1 / vb1 - 1 / vb2))/(1 / vb1 - 1 / vb2);
      end
      if ixe==1
       dc1 = (xsol (2) / 3  *(1 / vw - 1 / vb1) + (tetha / Vot))/ (1 / vw - 1 / vb1) ;  
       dc2 = dc1+(xsol (2)/3 * (1 / vb1) + (xsol (3) +xsol (4))/8*(1 / vb1 - 1 / (vb2)))/(1 / vb1 - 1 / vb2); 
      end
      
      if dc1<=dc2  
      % Case 1.1
       if svar<=dc1 % dc= dimension of walking
        kk=11;
        Acw = 2 * svar * svar;
        dcw = 2*(svar + svar) / 3;
        Acb1 = 0;
        dcb1 = 0;
        Acb2 = 0;
        dcb2 = 0;
       end 

       %Case 2
       if dc1 <= svar & dc2 >= svar  
        kk=12;
        Acw = 2 * dc1 * dc1;
        dcw = (dc1 + dc1) / 3;
        Acb1 = 2*svar*svar - 2*dc1 * dc1;
        dcb1 = 2*(2*svar*svar*svar/3-4*dc1*dc1*dc1/3)/Acb1;
        Acb2 = 0;
        dcb2 = 0;
       end 

       %Case 1.3
       if   dc2 <= svar  
        kk=13;
        Acw = 2 * dc1 * dc1;
        dcw = 2*(dc1 + dc1) / 3;
        Acb1 = 2*dc2*dc2-2*dc1*dc1;
        dcb1 =2/Acb1*(4*dc2*dc2*dc2/3-4*dc1*dc1*dc1/3);
        Acb2 = 2*svar * svar - 2*dc2*dc2;
        dcb2 =2*(4*svar*svar*svar/3-4*dc2*dc2*dc2/3)/Acb2; 
       end 
       
       A_cwT = A_cwT + (fc12(i, j) + fc3(i, j)) * Acw * dcw / vw / svar / svar/2* deltax * Deltay;
       A_cb1 = A_cb1 + (fc12(i, j) + fc3(i, j)) * Acb1 * dcb1 / vb1 / svar / svar/2 * deltax * Deltay;
       IVTT_c = IVTT_c + (fc12(i, j) + fc3(i, j)) * Acb2 * dcb2 / vb2 / svar / svar/2 * deltax * Deltay;
       q_c_uno = deltax*Deltay * Acb2 / svar / svar/2*(fc_primer (i, j) + fc_segundo (i, j));
       q_c_dos = deltax*Deltay * Acb2 / svar / svar/2*(fc_tercer (i, j) + fc_cuarto (i, j));

       if q_c_uno > q_c_dos 
          q_c = q_c_uno;
       else
          q_c = q_c_dos;
       end 


       M2 = 0;
       term1= (fb_out (i, j) + fc_primer (i, j) + fc_segundo (i, j));
       term2= (fb_in (i, j) + fc_tercer (i, j) + fc_cuarto (i, j));
       if term1>=term2  
         M2 = (term1-term2) * horizon * deltax * Deltay;
       end 

       Zoper=0;
       Zrepoteams=0; % Repoteams= Number of rebalancing teams throughout the city
       if term1>=term2 
          Linehaultrips= (term1-term2) * horizon * Dx * Dy*2/vancap;
          Distlinehaul=Linehaultrips*factordis*(Dx*Dy) ^ (1/2);
          Distpeddling=1.1*Dx*Dy*(1/xsol (2)/xsol (2)) ^ (1/2);
          Repoteams= ((Distlinehaul+Distpeddling)/vvan+2*deltaload*(Dx*Dy*(term1*horizon*2/xsol (2)/xsol (2)) ^ (1/2) +Dx*Dy*(term1-term2) * horizon))/horizon;        
          Repoteamscell=Repoteams/40/40; % Number of rebalancing teams divided by study cell area          
          Zrepoteams=Repoteamscell*Crepoteams; % Rebalancing cost depending on the number of teams
       end

       Zoper=term2*deltax*Deltay*Coper; % Operating cost of bike trips
       M3 =1.28* sqrt (2 * (Bpax+Cpax)/Dx/Dy*horizon* 2 / xsol (2) /xsol (2)) * deltax * Deltay;
       M4 =1.28* sqrt (IVTT_bx + IVTT_by + IVTT_c+A_cb1+A_bb;
       M1=IVTT_bx + IVTT_by + IVTT_c+A_cb1+A_bb; 

       % I consider stock in bike stations in the periphery
       Npt = 0;
       if i < nx_out | i > nx_in + nx_out | j < ny_out | j > ny_in + ny_out 
           Npt = deltax * Deltay * deltapt; 
       end 
       M = IVTT_bx + IVTT_by + IVTT_c + M2+ M3+ M4; % M= fleet size (number of vehicles)
       if fb (i, j) ==0&fc12(i, j) ==0&fc3(i, j) ==0
          Ns=0;
          Npt=0;
          xsol (2) =Dx;
       end
      
       % FORMULATE OBJECTIVE FUNCTION
       Z = Zoper+Zrepoteams+(Cm * M + Cl * L + Cs * (Npt + Ns) + Vot * (IVTT_ax + IVTT_ay + IVTT_bx + IVTT_by + IVTT_c + A_ab + A_bw + A_bb + A_cb1+A_cwT)); 
           
       kkk=kkk+1;
               
       % CAPACITY
       if i < nx_out | i > nx_in + nx_out | j < ny_out | j > ny_in + ny_out 
             qxpos = q_axpos (i, j) + q_bxpos (i, j) + q_c;
             qxneg = q_axneg (i, j) + q_bxneg (i, j) + q_c;
             qypos = q_aypos (i, j) + q_bypos (i, j) + q_c;
             qyneg = q_ayneg (i, j) + q_byneg (i, j) + q_c;
       else
             qxpos = q_axpos (i, j) + q_bxpos (i, j);
             qxneg = q_axneg (i, j) + q_bxneg (i, j);
             qypos = q_aypos (i, j) + q_bypos (i, j);
             qyneg = q_ayneg (i, j) + q_byneg (i, j);
       end 

       if qxpos<=C & qxneg<=C & qypos<=C& qyneg<=C &Z<=Zmin
           Zmin=Z;
           Zdef (i, j) =Z;
           ss (i, j) =xsol (2);
           sBx (i, j) =xsol (3);
           sBy (i, j) =xsol (4);
           Mdef (i, j) =M;
           Ldef (i, j) =L;
           Nsdef (i, j) =Ns+Npt;
           IVTT_adef (i, j) =IVTT_ax +IVTT_ay;
           IVTT_bdef (i, j) =IVTT_bx +IVTT_by;
           IVTT_cdef (i, j) =IVTT_c;
           A_abdef (i, j) =A_ab;
           A_bwdef (i, j) =A_bw;
           A_bbdef (i, j) =A_bb;
           A_cb1def (i, j) =A_cb1;
           A_cwTdef (i, j) =A_cwT;
           Zudef (i, j) =Vot * (IVTT_ax + IVTT_ay + IVTT_bx + IVTT_by + IVTT_c + A_ab + A_bw + A_bb + A_cb1+A_cwT);
           Zadef (i, j) = (Cm * M + Cl * L + Cs * (Npt + Ns));
           Zoperdef (i, j) =Zoper;
           Zrepodef (i, j) =Zrepoteams;
           Tdef (i, j) =IVTT_ax + IVTT_ay + IVTT_bx + IVTT_by + IVTT_c + A_ab + A_bw + A_bb + A_cb1+A_cwT;

           qxposdef (i, j) = qxpos;
           qxnegdef (i, j) = qxneg;  
           qyposdef (i, j) = qypos; 
           qynegdef (i, j) = qyneg;
       end
     end
    end
    end
  end  
  end
  end
          qxposdef
          qxnegdef  
          qyposdef 
          qynegdef
          
% DISCRETIZATION OF THE NETWORK IN CORRIDORS
 
% horizontal axes lower quadrant
 yd_inf=zeros (20,300);  
 ncorredorx_inf=zeros (300);
 ycorr_inf=zeros (20,300);
 nlinkfw_inf=zeros (20,300);
 nlinkbw_inf=zeros (20,300);

 for j=1:40
     yd_inf (j,1) =0;
     ncorredorx_inf (j) =1;
     disw=0; % disw= distance traveled
     dmin=0;
     for ii=1:20
         i=40-ii+1;
         ii;
         k=0;
         while k==0
             if disw+(Deltay*ii-dmin)/sBy (i, j) <1
               k=1;
             else
               k=0;
               ncorredorx_inf (j) =ncorredorx_inf (j) +1;
               yd_inf (j, ncorredorx_inf (j)) =dmin+sBy (i, j) *(1-disw);
               disw=0;
              
               dmin=yd_inf (j, ncorredorx_inf (j));
             end    
         end
         dmin=Deltay*ii;
         disw=min (disw+Deltay/sBy (i, j), Disw+(Deltay*ii-yd_inf (j, ncorredorx_inf (j)))/sBy (i, j));
     end

     % adjust of last lane distance 
     if Deltay*20-yd_inf (j, ncorredorx_inf (j))>=0.5*sBy (i, j)
         ncorredorx_inf (j) =ncorredorx_inf (j) +1;
         yd_inf (j, ncorredorx_inf (j)) =Deltay*20;
     else
         yd_inf (j, ncorredorx_inf (j)) =Deltay*20;
     end
 end

 % Calculate the y-coordinate of the corridor that passes through the middle of the coverage zone
 ncorxmax_inf=0; % max number of lanes in lower x quadrant
 for j=1:40
     if ncorredorx_inf (j)>ncorxmax_inf
         ncorxmax_inf=ncorredorx_inf (j);
     end
     for k=1: ncorredorx_inf (j)-1
         ycorr_inf (j, k) =yd_inf (j, k) + (yd_inf (j, k+1)-yd_inf (j, k))/2;
     end
  end

  % join rail sections with forks
 for j=1:39
     for k=1: ncorredorx_inf (j)
         dmax=10000000;
         for kk=1: ncorredorx_inf (j+1)
             if abs (ycorr_inf (j, k)-ycorr_inf (j+1, kk)) < dmax
                 nlinkfw_inf (j, k) =kk; % number of forward links lower quadrant
                 dmax=abs (ycorr_inf (j, k)-ycorr_inf (j+1, kk));
             end
         end
     end
 end
 
 for j=2:40
     for k=1: ncorredorx_inf (j)
         dmax=10000000;
         for kk=1: ncorredorx_inf (j-1)
             if abs (ycorr_inf (j, k)-ycorr_inf (j-1, kk)) < dmax
                 nlinkbw_inf (j, k) =kk; % number of backward links lower quadrant
                 dmax=abs (ycorr_inf (j, k)-ycorr_inf (j-1, kk));
             end
         end
     end
 end
  
 % horizontal axes upper quadrant
 yd_sup=zeros (20,300);
 yd_aux=zeros (20,300);
 ncorredorx_sup=zeros (300);
 ycorr_sup=zeros (20,300);
 nlinkfw_sup=zeros (20,300);
 nlinkbw_sup=zeros (20,300);
 
for j=1:40
     yd_sup (j,1) =0;
     ncorredorx_sup (j) =1;
     disw=0;
     dmin=0;
   
  for i=1:20
        k=0;
         while k==0
             if disw+ (Deltay*i-dmin)/sBy (i, j) <1
               k=1;
             else
               k=0;
               ncorredorx_sup (j) =ncorredorx_sup (j) +1;
               yd_sup (j, ncorredorx_sup (j)) =dmin+sBy (i, j) *(1-disw);
               disw=0;
               dmin=yd_sup (j, ncorredorx_sup (j));
             end    
         end
         dmin=Deltay*i;
         disw=min (disw+Deltay/sBy (i, j), Disw+ (Deltay*i-yd_sup (j, ncorredorx_sup (j)))/sBy (i, j));
     end

     % adjust of last lane distance
     if Deltay*20-yd_sup (j, ncorredorx_sup (j))>=0.5*sBy (i, j)
         ncorredorx_sup (j) =ncorredorx_sup (j) +1;
         yd_sup (j, ncorredorx_sup (j)) =Deltay*20;
     else
         yd_sup (j, ncorredorx_sup (j)) =Deltay*20;
     end
 end
 
for j=1:40
     for kkk=1: ncorredorx_sup (j)
         yd_aux (j, kkk) =yd_sup (j, kkk);
         yd_sup (j, kkk) =0;
     end
     for kkk=1: ncorredorx_sup (j)
         yd_sup (j, kkk) =Deltay*40-yd_aux (j, kkk);
     end
 end
 
 % Calculate the y-coordinate of the corridor that passes through the middle of the coverage zone
 ncorxmax=0;

 for j=1:40
     if ncorredorx_inf (j) +ncorredorx_sup (j)>ncorxmax
         ncorxmax=ncorredorx_inf (j) +ncorredorx_sup (j);
     end
     for k=1: ncorredorx_inf (j)-1
         ycorr_inf (j, k) =yd_inf (j, k) + (yd_inf (j, k+1)-yd_inf (j, k))/2;
     end
     
for k=1: ncorredorx_sup (j)-1
         ycorr_sup (j, k) =yd_sup (j, k) -(-yd_sup (j, k+1) +yd_sup (j, k))/2;
     end
 end

 ycorr_inf
 ycorr_sup
 ncorxmax

 % join rail sections with forks
 for j=1:39
     for k=1: ncorredorx_inf (j)
         dmax=10000000;
         for kk=1: ncorredorx_inf (j+1)
             if abs (ycorr_inf (j, k)-ycorr_inf (j+1, kk)) < dmax
                 nlinkfw_inf (j, k) =kk;
                 dmax=abs (ycorr_inf (j, k)-ycorr_inf (j+1, kk));
             end
         end
     end
 end
 
 for j=1:39
     for k=1: ncorredorx_sup (j)
         dmax=10000000;
         for kk=1: ncorredorx_sup (j+1)
             if abs (ycorr_sup (j, k)-ycorr_sup (j+1, kk)) < dmax
                 nlinkfw_sup (j, k) =kk;
                 dmax=abs (ycorr_sup (j, k)-ycorr_sup (j+1, kk));
             end
         end
     end
 end

 for j=2:40
     for k=1: ncorredorx_inf (j)
         dmax=10000000;
         for kk=1: ncorredorx_inf (j-1)
             if abs (ycorr_inf (j, k)-ycorr_inf (j-1, kk)) < dmax
                 nlinkbw_inf (j, k) =kk;
                 dmax=abs (ycorr_inf (j, k)-ycorr_inf (j-1, kk));
             end
         end
     end
 end

 for j=2:40
     for k=1: ncorredorx_sup (j)
         dmax=10000000;
         for kk=1: ncorredorx_sup (j-1)
             if abs (ycorr_sup (j, k)-ycorr_sup (j-1, kk)) < dmax
                 nlinkbw_sup (j, k) =kk;
                 dmax=abs (ycorr_sup (j, k)-ycorr_sup (j-1, kk)); % y coordinate of point k of corridor j in the upper quadrant
             end
         end
     end
 end
 
 figure 
 hold on % To overlay multiple graphs on the same figure

 % Define boundaries and aspect ratio
 xlim ([0, Dx]);
 ylim ([0, Dy]);
 daspect ([1,1,1]);

 % Draw lines and corridors
 plot ([0, Dx], [Dy, Dy],'k') % Draw a horizontal line at the top of the quadrant
 plot ([Dx, Dx], [0, Dy],'k') % Draw a vertical line on the right side of the quadrant
 
 for j=1:40
  for k=1: ncorredorx_inf (j)
      plot ([deltax*(j-1), Deltax*(j)], [ycorr_inf (j, k),ycorr_inf (j, k)],'b')
  end
  for k=1: ncorredorx_sup (j)
      plot ([deltax*(j-1), Deltax*(j)], [ycorr_sup (j, k),ycorr_sup (j, k)],'b')
  end
 end  
 for j=2:40
    for k=1: ncorredorx_inf (j)
     plot ([deltax*(j-1), Deltax*(j-1)], [ycorr_inf (j-1,nlinkbw_inf (j, k)),ycorr_inf (j, k)],'c')
    end 
    for k=1: ncorredorx_sup (j)
     plot ([deltax*(j-1), Deltax*(j-1)], [ycorr_sup (j-1,nlinkbw_sup (j, k)),ycorr_sup (j, k)],'c')
    end 
 end
 for j=1:39
    for k=1: ncorredorx_inf (j)
     plot ([deltax*(j), Deltax*(j)], [ycorr_inf (j+1,nlinkfw_inf (j, k)),ycorr_inf (j, k)],'c')
    end 
   
 for k=1: ncorredorx_sup (j)
     plot ([deltax*(j), Deltax*(j)], [ycorr_sup (j+1,nlinkfw_sup (j, k)),ycorr_sup (j, k)],'c')
    end 
 end
 hold off
    
  % LEFT VERTICAL AXES 
 xd_left=zeros (20,300);
 ncorredory_left=zeros (300);
 xcorr_left=zeros (20,300);
 xnlinkfw_left=zeros (20,300);
 xnlinkbw_left=zeros (20,300);
 for i=1:40
     xd_left (i,1) =0;
     ncorredory_left (i) =1;
     disw=0;
     dmin=0;
     for j=1:20
         k=0;
         while k==0
             if disw+ (deltax*j-dmin)/sBx (i, j) <1
               k=1;
             else
               k=0;
               ncorredory_left (i) =ncorredory_left (i) +1;
               
               xd_left (i, ncorredory_left (i)) =dmin+sBx (i, j) *(1-disw);
               disw=0;              
               dmin=xd_left (i, ncorredory_left (i));
             end    
         end
         dmin=deltax*j;
         disw=min (disw+deltax/sBx (i, j), Disw+ (deltax*j-xd_left (i, ncorredory_left (i)))/sBx (i, j));
     end
  
   % adjust of last lane distance
     if deltax*20-xd_left (i, ncorredory_left (i))>=0.5*sBx (i, j)
         ncorredory_left (i) =ncorredory_left (i) +1;
         xd_left (i, ncorredory_left (i)) =deltax*20;
     else
         xd_left (i, ncorredory_left (i)) =deltax*20;
     end
 end

 % Calculate the x-coordinate of the lane that passes through the middle of the coverage zone
 ncorymax_left=0;
 for i=1:40
     if ncorredory_left (i)>ncorymax_left
         ncorymax_left=ncorredory_left (i);
     end
     for k=1: ncorredory_left (i)-1
         xcorr_left (i, k) =xd_left (i, k) + (xd_left (i, k+1)-xd_left (i, k))/2;
     end
  end

 % join rail sections with forks
 for i=1:39
     for k=1: ncorredory_left (i)
         dmax=10000000;
         for kk=1: ncorredory_left (i+1)
             if abs (xcorr_left (i, k)-xcorr_left (i+1, kk)) < dmax
                 xnlinkfw_left (i, k) =kk;
                 dmax=abs (xcorr_left (i, k)-xcorr_left (i+1, kk));
             end
         end
     end
 end

 for i=2:40
     for k=1: ncorredory_left (i)
         dmax=10000000;
         for kk=1: ncorredory_left (i-1)
             if abs (xcorr_left (i, k)-xcorr_left (i-1, kk)) < dmax
                xnlinkbw_left (i, k) =kk;
                dmax=abs (xcorr_left (i, k)-xcorr_left (i-1, kk));
             end
         end
     end
 end

 % RIGHT VERTICAL AXES
 xd_rig=zeros (20,300);
 xd_aux=zeros (20,300);
 ncorredory_rig=zeros (300);
 xcorr_rig=zeros (20,300);
 xnlinkfw_rig=zeros (20,300);
 xnlinkbw_rig=zeros (20,300);

 for i=1:40
     xd_rig (i,1) =0;
     ncorredory_rig (i) =1;
     disw=0;
     dmin=0;
     for jj=1:20
         j=40-jj+1;
         k=0;
         while k==0
             if disw+ (deltax*jj-dmin)/sBx (i, j) <1
               k=1;
             else
               k=0;
               ncorredory_rig (i) =ncorredory_rig (i) +1;
               xd_rig (i, ncorredory_rig (i)) =dmin+sBx (i, j) *(1-disw);
               disw=0;
               dmin=xd_rig (i, ncorredory_rig (i));
             end    
         end
         dmin=deltax*jj;
         disw=min (disw+deltax/sBx (i, j), Disw+ (deltax*jj-xd_rig (i, ncorredory_rig (i)))/sBx (i, j));
     end     
   % adjust of last lane distance 
     if deltax*20-xd_rig (i, ncorredory_rig (i))>=0.5*sBx (i, j)
         ncorredory_rig (i) =ncorredory_rig (i) +1;
         xd_rig (i, ncorredory_rig (i)) =deltax*20;
     else
         xd_rig (i, ncorredory_rig (i)) =deltax*20;
     end
 end

 for i=1:40
     for kkk=1: ncorredory_rig (i)
         xd_aux (i, kkk) =xd_rig (i, kkk);
         xd_rig (i, kkk) =0;
     end
     for kkk=1: ncorredory_rig (i)
         xd_rig (i, kkk) =deltax*40-xd_aux (i, kkk);
     end
 end

% I calculate the x-coordinate of the lane that passes through the middle of the coverage zone
 ncorymax=0;
 for i=1:40
     if ncorredory_left (i) +ncorredory_rig (i)>ncorymax
         ncorymax=ncorredory_left (i) +ncorredory_rig (i);
     end
     for k=1: ncorredory_left (i)-1
         xcorr_left (i, k) =xd_left (i, k) + (xd_left (i, k+1)-xd_left (i, k))/2;
     end
     for k=1: ncorredory_rig (i)-1
         xcorr_rig (i, k) =xd_rig (i, k) -(-xd_rig (i, k+1) +xd_rig (i, k))/2;
     end
 end

% join rail sections with forks
 for i=1:39
     for k=1: ncorredory_left (i)
         dmax=10000000;
         for kk=1: ncorredory_left (i+1)
             if abs (xcorr_left (i, k)-xcorr_left (i+1, kk)) < dmax
                 xnlinkfw_left (i, k) =kk;
                 dmax=abs (xcorr_left (i, k)-xcorr_left (i+1, kk));
             end
         end
     end

     for k=1: ncorredory_rig (i)
         dmax=10000000;
         for kk=1: ncorredory_rig (i+1)
             if abs (xcorr_rig (i, k)-xcorr_rig (i+1, kk)) < dmax
                 xnlinkfw_rig (i, k) =kk;
                 dmax=abs (xcorr_rig (i, k)-xcorr_rig (i+1, kk));
             end
         end
     end
 end

 for i=2:40
     for k=1: ncorredory_left (i)
         dmax=10000000;
         for kk=1: ncorredory_left (i-1)
             if abs (xcorr_left (i, k)-xcorr_left (i-1, kk)) < dmax
                 xnlinkbw_left (i, k) =kk;
                 dmax=abs (xcorr_left (i, k)-xcorr_left (i-1, kk));
             end
         end
     end

     for k=1: ncorredory_rig (i)
         dmax=10000000;
         for kk=1: ncorredory_rig (i-1)
             if abs (xcorr_rig (i, k)-xcorr_rig (i-1, kk)) < dmax
                 xnlinkbw_rig (i, k) =kk;
                 dmax=abs (xcorr_rig (i, k)-xcorr_rig (i-1, kk));
             end
         end
     end
 end
 xcorr_left
 xcorr_rig
 ncorymax
 


 figure 
 hold on
 xlim ([0, Dx]);
 ylim ([0, Dy]);
 daspect ([1,1,1]);
 plot ([0, Dx], [Dy, Dy],'k')
 plot ([Dx, Dx], [0, Dy],'k')

 for i=1:40
  for k=1: ncorredory_left (i)
      plot ([xcorr_left (i, k), xcorr_left (i, k)], [Deltay*(i-1), Deltay*(i)],'r')
  end
  for k=1: ncorredory_rig (i)
      plot ([xcorr_rig (i, k), xcorr_rig (i, k)], [Deltay*(i-1), Deltay*(i)],'r')
  end
 end  

 for i=2:40
    for k=1: ncorredory_left (i)
     plot ([xcorr_left (i-1, xnlinkbw_left (i, k)), xcorr_left (i, k)], [Deltay*(i-1), Deltay*(i-1)],'m')
    end 
    for k=1: ncorredory_rig (i)
     plot ([xcorr_rig (i-1, xnlinkbw_rig (i, k)), xcorr_rig (i, k)], [Deltay*(i-1), Deltay*(i-1)],'m')
    end 
 end

 for i=1:39
    for k=1: ncorredory_left (i)
     plot ([xcorr_left (i+1, xnlinkfw_left (i, k)), xcorr_left (i, k)], [Deltay*(i), Deltay*(i)],'m')
    end 
    for k=1: ncorredory_rig (i)
     plot ([xcorr_rig (i+1, xnlinkfw_rig (i, k)), xcorr_rig (i, k)], [Deltay*(i), Deltay*(i)],'m')
    end 
 end
 hold off


 % GENERATION OF THE JOINT GRAPH
 figure 
 hold on
 xlim ([0, Dx]);
 ylim ([0, Dy]);
 daspect ([1,1,1]);
 plot ([0, Dx], [Dy, Dy],'k')
 plot ([Dx, Dx], [0, Dy],'k')

 for i=1:40
  for k=1: ncorredory_left (i)
      plot ([xcorr_left (i, k), xcorr_left (i, k)], [Deltay*(i-1), Deltay*(i)],'r')
  end
    for k=1: ncorredory_rig (i)
      plot ([xcorr_rig (i, k), xcorr_rig (i, k)], [Deltay*(i-1), Deltay*(i)],'r')
    end
 end  

 for i=2:40
    for k=1: ncorredory_left (i)
     plot ([xcorr_left (i-1, xnlinkbw_left (i, k)), xcorr_left (i, k)], [Deltay*(i-1), Deltay*(i-1)],'m')
    end 
    for k=1: ncorredory_rig (i)
     plot ([xcorr_rig (i-1, xnlinkbw_rig (i, k)), xcorr_rig (i, k)], [Deltay*(i-1), Deltay*(i-1)],'m')
    end 
 end

 for i=1:39
    for k=1: ncorredory_left (i)
     plot ([xcorr_left (i+1, xnlinkfw_left (i, k)), xcorr_left (i, k)], [Deltay*(i), Deltay*(i)],'m')
    end 
    for k=1: ncorredory_rig (i)
     plot ([xcorr_rig (i+1, xnlinkfw_rig (i, k)), xcorr_rig (i, k)], [Deltay*(i), Deltay*(i)],'m')
    end 
 end


 for j=1:40
  for k=1: ncorredorx_inf (j)
      plot ([deltax*(j-1), Deltax*(j)], [ycorr_inf (j, k),ycorr_inf (j, k)],'b')
  end
  for k=1: ncorredorx_sup (j)
      plot ([deltax*(j-1), Deltax*(j)], [ycorr_sup (j, k),ycorr_sup (j, k)],'b')
  end
 end  
 for j=2:40
    for k=1: ncorredorx_inf (j)
     plot ([deltax*(j-1), Deltax*(j-1)], [ycorr_inf (j-1,nlinkbw_inf (j, k)),ycorr_inf (j, k)],'c')
    end 
    for k=1: ncorredorx_sup (j)
     plot ([deltax*(j-1), Deltax*(j-1)], [ycorr_sup (j-1,nlinkbw_sup (j, k)),ycorr_sup (j, k)],'c')
    end 
 end
 for j=1:39
    for k=1: ncorredorx_sup (j)
     plot ([deltax*(j), Deltax*(j)], [ycorr_sup (j+1,nlinkfw_sup (j, k)),ycorr_sup (j, k)],'c')
    end 
 end
 hold off
