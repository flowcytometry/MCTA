! 
!   MCTA v1.1
!   Copyright by Carine Beatrici, Fabrício A. B. Silva
!
!
! analised_colors......the colors that are going to be taken into account in the tendency 
! use_in_average..=..0 that color is not used in the average calculation
!               .....1 the i color is used to calculate the color tendency
!
!
!
!
program MCTA
implicit none

INTERFACE
real FUNCTION Median(X, N)
      REAL, DIMENSION(1:N), INTENT(IN)    :: X
      INTEGER, INTENT(IN)                :: N
     END FUNCTION Median     
END INTERFACE

! Variables declaration *************************************************************************************************
! File names and temporary files
character:: file_name*40,file_name_dat*40,file_name_txt*40
character:: anything*20,command*100,line*1000
integer:: rad_name,number_rows_original_file=0

! Events and counter variables
integer:: control=0,counter,event_accepted
integer:: i,j,k,vhmax,vhmin, number_values_median
integer:: total_events_accepted=0, total_events_rejected=0

! file rows informations 
integer:: rows,fluorescent_channels,color_begin,number_color_analysed
logical:: file_fcs_exist

! Calculated Values
real:: average_l,sum_x,sum_y,hue_med,light_med,x_color,somatorio,sumsqr
real:: total_light,max_light

real,allocatable:: values(:,:),average(:),in_data(:), median_values(:)
integer,allocatable:: lambda(:),analised_colors(:),use_in_average(:),choosen_rows(:),order(:)
real,allocatable:: hue(:),hue_real(:),light(:),compensation(:,:),mean(:),std(:),variance(:),median_result(:)

! Deleting old output files **********************************************************************************************
file_name=" "
file_name(1:8)="nofl.dat"
inquire( file=trim(file_name), exist=file_fcs_exist )
if(file_fcs_exist)call system("rm nofl.dat")
file_name=" "
file_name(1:8)="fluo.dat"
inquire( file=trim(file_name), exist=file_fcs_exist )
if(file_fcs_exist)call system("rm fluo.dat")
file_name=" "
file_name(1:16)="result_file.eps"
inquire( file=trim(file_name), exist=file_fcs_exist )
if(file_fcs_exist)call system("rm result_file.eps")
file_name=" "
file_name(1:19)="result_file-deb.eps"
inquire( file=trim(file_name), exist=file_fcs_exist )
if(file_fcs_exist)call system("rm result_file-deb.eps")

! Reading the parameters file *******************************************************************************************
open(20,file="in.dat")
file_name=" "
command=" "
read(20,*,ERR=1000,END=1000)file_name
inquire( file=trim(file_name), exist=file_fcs_exist )
if(file_fcs_exist)then
  print *, ''//achar(27)//"[32m FCS File found: ",file_name," "//achar(27)//"[0m"
else
  print *, ''//achar(27)//"[31m FCS File ",file_name," not found!!!"//achar(27)//"[0m"
  print *, "Please make sure the file",file_name," is in the current folder."
  stop
end if

file_name_txt=file_name
file_name_dat=file_name
do i=1,40
  if(file_name(i:i)=='.')then
    rad_name=i-1
    file_name(i:i+4)  =''
    file_name_txt(i:i+4)  ='.txt'
    file_name_dat(i:i+4)  ='.dat'
    exit
  end if
end do

goto 1001
1000 print *, ''//achar(27)//"[31m Error in file in.dat "//achar(27)//"[0m"
print *, "The file in.dat is required, to proceed create or edit the file "
print *, "Please include the name of the FCS file to be analised in the first line of the file in.dat"
stop
1001 continue

! converting the fcs file to a text temporary file **********************************************************************
command = "Rscript "//"FCS2CSV.R "//file_name;
call system(command) 
print *, ''//achar(27)//"[32m FCS File ",file_name(1:rad_name),".fcs successfully converted"//achar(27)//"[0m"

open(30,file=file_name_txt)
open(40,file="result_file")
line=''
read(30,'(a1000)',END=100,ERR=200)line
write(*,*)"Temporary text file: ",file_name_txt," "
write(*,'((a12)$)')"First line: "
do i=1,1000
  if(line(i:i).ne." ")write(*,'((a1)$)')line(i:i)
end do
write(*,*)' '


read(20,*,ERR=1002,END=1002)rows
goto 1003
1002 print *, ''//achar(27)//"[31m Error in file in.dat"//achar(27)//"[0m"
print *, " in.dat file incomplete "
print *, "Please inform the number of needed columns of the .fcs file"
stop
1003 continue

allocate(choosen_rows(rows))
choosen_rows=0
read(20,*,ERR=1004,END=1004)fluorescent_channels
read(20,*,ERR=1004,END=1004)choosen_rows
write(*,*)"Columns used:",choosen_rows
goto 1005
1004 print *, ''//achar(27)//"[31m Error in  file in.dat"//achar(27)//"[0m"
print *, "in.dat file incomplete "
print *, "Please check the .fcs file"
stop
1005 continue

!converting the file -  required data only ************************************************************
do i =1,1000
  if(line(i:i)==',')then
    line(i:i)=' '
    number_rows_original_file=number_rows_original_file+1
  end if
  if(line(i:i)=='"')line(i:i)=' '
end do
write(*,*)"Number of columns present in the input data file:",number_rows_original_file," "
allocate(in_data(number_rows_original_file))
!write(*,*)"first line of the file:",line
do while(control==0)
  line=''
  read(30,'(a1000)',END=100,ERR=200)line
  do i=1,1000
    if(line(i:i)==','.or.line(i:i)=='"')line(i:i)=' '
  end do
  write(40,*)line
end do 
200 write(*,*)"................................................................." 
print *, ''//achar(27)//"[31m There was an error in the file ",file_name_txt,".fcs"//achar(27)//"[0m"
print *, ''//achar(27)//"[31m Please check your fcs file and dependencies"//achar(27)//"[0m"
stop
100 print *, ''//achar(27)//"[32m FCS file processed"//achar(27)//"[0m"
command="cp result_file "//file_name_txt
call system(command)
close(40)
call system("rm result_file")
close(30)

!reopening the files ***************************************************************************************************
open(30,file=file_name_txt)
open(40,file=file_name_dat)
do while(control==0)
  line=''
  read(30,*,END=400,ERR=300)in_data
  do i =1,rows
    write(40,"(EN15.6E2,3x)",advance="no")in_data(choosen_rows(i))
  end do
  write(40,*)" "
end do
300 print *, ''//achar(27)//"[31m There was an error in the ",file_name_txt,".fcs"//achar(27)//"[0m"
print *, ''//achar(27)//"[31m Please check your fcs file and dependencies"//achar(27)//"[0m"


stop
400 print *, ''//achar(27)//"[32m FCS file processed"//achar(27)//"[0m"
write(*,*)"................................................................."
write(*,*)"................................................................." 
close(30)
close(40)

! allocating analisys data ********************************************************************************************
allocate(lambda(fluorescent_channels),hue(fluorescent_channels),light(fluorescent_channels),hue_real(fluorescent_channels))
allocate(mean(fluorescent_channels),std(fluorescent_channels),variance(fluorescent_channels),median_result(fluorescent_channels))
allocate(order(fluorescent_channels))
allocate(compensation(fluorescent_channels,fluorescent_channels))
read(20,*,ERR=1004,END=1004)lambda
allocate(average(fluorescent_channels),analised_colors(fluorescent_channels))
allocate(use_in_average(fluorescent_channels))
read(20,*,ERR=1004,END=1004)average
do i=1,fluorescent_channels
  read(20,*,ERR=1004,END=1004)compensation(i,1:fluorescent_channels)
  compensation(i,i)=0.0
end do
analised_colors=0
read(20,*,ERR=1004)analised_colors
close(20)

! Final in.dat reading ************************************************************************************************ 

color_begin = rows - fluorescent_channels
write(*,'((a8),10(i3,3x))')"lambda: ",lambda
do i=1,fluorescent_channels
   write(*,'((a2),10(f6.2,3x))')"C:",compensation(i,:)
end do
compensation = compensation/100.00
write(*,*)"................................................................." 
write(*,*)"Number of columns used: ",rows," " 
write(*,*)"Number of fluorescent channels chosen: ",fluorescent_channels," " 
!write(*,*)"what row in the file the color begins: ",color_begin,".................................." 
write(*,'((a40),10(f10.2,3x)$)')"Negative background cut-off values: ",average
write(*,*)"Input fcs data file: ",file_name_dat," " 

! counting the total number of events *************************************************************************************
write(*,*)"Counting the number of events" 
open(10,file=file_name_dat)
counter=0
do while(control==0)
  read(10,*,END=20)anything
  counter=counter+1
end do
20 continue
if (counter == 0) then
  print *, ''//achar(27)//"[31m ERROR: There are no events in the file",file_name,".fcs !!!!!!!!!"//achar(27)//"[0m"
  print *, ''//achar(27)//"[31m Please check the script FCS2CSV.R and the FCS file"//achar(27)//"[0m"
  stop
end if
if (counter<=20000) then
        print *, ''//achar(27)//"[32m Number of events in the FCS file .......",counter," "//achar(27)//"[0m"
else
  if (counter <=100000)then
    print *, ''//achar(27)//"[33m Number of events in the FCS file .......",counter," "//achar(27)//"[0m"
  else
    if (counter <= 1000000)then
      print *, ''//achar(27)//"[33m Number of events in the FCS file .......",counter," "//achar(27)//"[0m"
    else
      print *, ''//achar(27)//"[33m WARNING: Number of events in the FCS file .......",counter," "//achar(27)//"[0m"
    end if
  end if
end if

close(10)
open(10,file=file_name_dat)

! Memory data alocation ***************************************************************************************************
allocate(values(counter,rows))
do i=1,counter
  read(10,*)values(i,:)
end do
close(10)
! Applying data compensation *********************************************************************************************
write(*,*)"Fluorescence data compensation"
do i=1,counter
  do j=1,fluorescent_channels
    do k=1,fluorescent_channels
      values(i,color_begin+j)=values(i,color_begin+j)-compensation(k,j)*values(i,color_begin+k)
    end do
  end do
end do
write(*,*)"................................................................." 
! Removing the background values of the data *******************************************************************************
write(*,*)"Removing background cut-off values:",average," " 
do i=1,counter
  do j=1,fluorescent_channels
    values(i,j+color_begin)=values(i,j+color_begin)-average(j)
    if(values(i,j+color_begin)<0.0)then
      values(i,j+color_begin)=0.0
    end if
  end do
end do
where(values<0.0)values=0.0
write(*,*)"................................................................." 
! Saving the compensated data for other analisys ***************************************************************************
file_name(rad_name+1:rad_name+8)='-cmp.dat'
write(*,*)"Saving compensated data: ",file_name," " 
open(10,file=file_name)
do i=1,counter
  write(10,*)values(i,:)
end do
close(10)
write(*,*)"................................................................." 

! Building the base vectors
write(*,*)"Building base vectors" 
vhmax = 0
vhmin = 100000
do j=1,fluorescent_channels
    if (lambda(j)>vhmax) vhmax=lambda(j)
    if (lambda(j)<vhmin) vhmin=lambda(j)
end do
! to avoid any color to reach the end of the scale
vhmax = vhmax + 0.02*vhmax
vhmin = vhmin - 0.02*vhmin
! calculation adjusted for the color bar of gnuplot -> (vhmax - lambda(j))
do j=1,fluorescent_channels
  hue_real(j)=((vhmax - lambda(j))*(6.2830))/(vhmax-vhmin)
end do
write(*,*)"hue real (radians):",hue_real
write(*,*)"hue real (degrees):",hue_real/6.2830*360
write(*,*)"................................................................."
!Sorting
do j=1,fluorescent_channels
    order(j) = j
    hue(j) = hue_real(j)
end do
!Insertion sort
do i = 2, fluorescent_channels
        x_color = hue(i)
        !k = ordem(i)
        j = i - 1
        do while (j >= 1)
            if (hue(j) <= x_color) exit
            hue(j + 1) = hue(j)
            k = order(j)
            order(j) = order(j+1)
            order(j+1) = k
            j = j - 1
        end do
        hue(j + 1) = x_color
        !order(i)=j+1
end do
!write(*,*)"order:",order
!write(*,*)"hue:",hue
! Distribute channels evenly along the hue circle
do j=1,fluorescent_channels
    hue(order(j)) = j*(6.2830/fluorescent_channels)
end do
!write(*,*)"order:",order
write(*,*)"hue recalculated (radians):",hue
write(*,*)"hue recalculated (degrees):",hue/6.2830*360
write(*,*)"................................................................."


! The output files
open(30,file="nofl.dat")
open(20,file="fluo.dat")


number_color_analysed = 0
do i=1,fluorescent_channels
   if(analised_colors(i)/=0)number_color_analysed=number_color_analysed+1
end do
write(*,*)"Number of colors excluded in the color tendency calculation:", number_color_analysed
write(*,*)"List of colors removed from tendency calculation:", analised_colors
use_in_average=1
do i=1,fluorescent_channels
  if(analised_colors(i)/=0)use_in_average(abs(analised_colors(i)))=0
end do
!write(*,*)"vetor use_in_average:", use_in_average

! Colour calculation for each event ********************************************************************************
write(*,*)"Please wait for the calculation"
write(*,*)"For larger files this can take some time"
do i=1,counter
  if(mod(i,int(counter/10))==0)print *, ''//achar(27)//'[32m',int(i/(counter/10))*10,'% '//achar(27)//'[0m'!write(*,*)int(i/(counter/10))*10,"% ....................."
  ! Calculo do L de cada evento e de cada canal de flourescencia
  if(number_color_analysed==0.or.number_color_analysed==fluorescent_channels)then
    ! all channels considered 
    average_l = sum(values(i,color_begin:rows))/real(fluorescent_channels) 
    do j=1,fluorescent_channels
      light(j) = values(i,j+color_begin)
    end do
    ! vectorial sum
    sum_x=0.0
    sum_y=0.0
    do j=1,fluorescent_channels
      sum_x=sum_x+light(j)*cos(hue(j))
      sum_y=sum_y+light(j)*sin(hue(j))
    end do
    light_med = sqrt(sum_x**2+sum_y**2)
    hue_med = atan2(sum_y,sum_x)
    !output
    do while(hue_med<0)
       hue_med=hue_med+6.2830
    end do
    do while(hue_med>=6.2830)
       hue_med=hue_med-6.2830
    end do
    hue_med = (hue_med/6.2830) * 360.0
    max_light=0.0
    total_light=0.0
    do j=1,fluorescent_channels
      if(values(i,color_begin+j)>max_light)max_light=values(i,color_begin+j)
      total_light=total_light+values(i,color_begin+j)
    end do

    if(light_med<=1.0d-5.or.isnan(hue_med))then
      write(30,*)values(i,1:color_begin),hue_med,1.0,0.5,values(i,color_begin:rows)
      total_events_rejected = total_events_rejected + 1
    else
      write(20,*)values(i,1:color_begin),hue_med,1.0,0.5,max_light/total_light,values(i,color_begin:rows)
      total_events_accepted = total_events_accepted + 1
    end if
  else
    ! one or more fluorescent_channels beeing analised
    average_l=0.0
    do j=1,fluorescent_channels
       average_l = average_l + values(i,color_begin+j)*use_in_average(j)
    end do
    average_l = average_l/real(fluorescent_channels-number_color_analysed+1) 
    do j=1,fluorescent_channels
       light(j) = values(i,j+color_begin)*use_in_average(j)
    end do
    sum_x=0.0
    sum_y=0.0
    do j=1,fluorescent_channels
      sum_x=sum_x+light(j)*cos(hue(j))
      sum_y=sum_y+light(j)*sin(hue(j))
    end do
    hue_med = atan2(sum_y,sum_x)
    light_med = sqrt(sum_x**2+sum_y**2)
    do while(hue_med < 0.0)
      hue_med = hue_med + 6.2830
    end do
    do while(hue_med > 6.2830)
      hue_med = hue_med - 6.2830
    end do
    hue_med = (hue_med/6.2830) * 360
    ! the new test: if the analised_colors is positive test for positive if the analised_colors is negative test for negative 
    event_accepted=1
    do k=1,fluorescent_channels
      ! event_accepted - AND operator
      if(analised_colors(k)>0)then
        if(values(i,color_begin+analised_colors(k))<=0.0)then
          event_accepted=0
        end if
      else
        if(analised_colors(k)<0)then
          if(values(i,color_begin+abs(analised_colors(k)))>0.0)then
            event_accepted=0
          end if
        end if
      end if
    end do
    max_light=0.0
    total_light=0.0
    do j=1,fluorescent_channels
      if((values(i,color_begin+j)*use_in_average(j))>max_light)max_light=(values(i,color_begin+j)*use_in_average(j))
      total_light=total_light+(values(i,color_begin+j)*use_in_average(j))
    end do
    if(total_light<=0.0)total_light=1.0 
    if(event_accepted==0.or.isnan(hue_med).or.light_med<=1.0d-5)then
      write(30,*)values(i,1:color_begin),hue_med,1.0,0.5,values(i,color_begin:rows)
      !write(*,*)values(1:color_begin)
    else
      write(20,*)values(i,1:color_begin),hue_med,1.0,0.5,max_light/total_light,values(i,color_begin:rows)
      !write(*,*)values(1:color_begin)
    end if
  end if
end do

! compute channel statistics
allocate(median_values(counter))
do j=1,fluorescent_channels
  somatorio = 0.0
  sumsqr = 0.0
  number_values_median = 0
  if (use_in_average(j)==1)then  
    do i=1,counter
      event_accepted=1
      do k=1,fluorescent_channels
      ! event_accepted - AND operator
        if(analised_colors(k)>0)then
           if(values(i,color_begin+analised_colors(k))<=0.0)then
             event_accepted=0
           end if
        else
          if(analised_colors(k)<0)then
            if(values(i,color_begin+abs(analised_colors(k)))>0.0)then
              event_accepted=0
            end if
          end if
        end if
      end do

      if (event_accepted==1) then
              if (values(i,j+color_begin)/=0.0) then 
                 somatorio = somatorio + LOG10(values(i,j+color_begin))
                 sumsqr = sumsqr + LOG10(values(i,j+color_begin))*LOG10(values(i,j+color_begin))
                 number_values_median = number_values_median + 1
                 median_values(number_values_median) = values(i,j+color_begin)
              end if
      end if
    end do
    mean(j) = 10**(somatorio/number_values_median)
    variance(j) = (sumsqr - somatorio*somatorio/number_values_median)/(number_values_median-1)
    std(j) = 10**SQRT(variance(j))
    !write(*,*)"XXKKKK: ",median_values(1)
    median_result(j) = Median(median_values,number_values_median)
    write(*,*)"Channel:",j," Geometric Mean: ",mean(j),"Geo Std dev: ",std(j), "Median: ",median_result(j)
  else
    mean(j) = 0
    variance(j) = 0
    std(j) = 0
  end if
end do


! deleting temporary files
write(*,*)"................................................................."
write(*,*)"Deleting temporary files"
command = "rm "//file_name_txt
call system(command)
command = "rm "//file_name_dat
call system(command)
write(*,*)"................................................................."
write(*,*)"Generating images"
!write(*,*)" Press any key to continue"
call system("gnuplot script-gnu-color")
call system("gnuplot script-gnu-bw")
deallocate(values)
write(*,*)"................................................................."
write(*,*)"Execution finished successfully!"
!print *, ''//achar(27)//"[95m .......................THE......END.............................. "//achar(27)//"[0m"
!print *, ''//achar(27)//"[95m ...........Execution finished with success!!!.................... :)"//achar(27)//"[0m"

end program


! --------------------------------------------------------------------
! REAL FUNCTION  Median() :
!    This function receives an array X of N entries, copies its value
! to a local array Temp(), sorts Temp() and returns the median.
! --------------------------------------------------------------------

real function  Median(X, N)
      IMPLICIT  NONE
      INTEGER, INTENT(IN)                :: N
      REAL, DIMENSION(1:N), INTENT(IN)   :: X
      REAL, DIMENSION(1:N)               :: Temp
      INTEGER                            :: i,j
      REAL                               :: xi


      DO i = 1, N                       ! make a copy
         Temp(i) = X(i)
      END DO

      CALL  Sort(Temp, N)               ! sort the copy
      IF (MOD(N,2) == 0) THEN           ! compute the median
         Median = (Temp(N/2) + Temp(N/2+1)) / 2.0
      ELSE
         Median = Temp(N/2+1)
      END IF


   END function  Median

!--------------------------------------------------------------------
!SUBROUTINE Sort()
!   A Simple sorting subroutine
!   Implements the INSERTION SORT algorithm
!--------------------------------------------------------------------

   subroutine Sort(a, n)
    implicit none
    integer :: n, i, j
    real :: a(n), x

    do i = 2, n
        x = a(i)
        j = i - 1
        do while (j >= 1)
            if (a(j) <= x) exit
            a(j + 1) = a(j)
            j = j - 1
        end do
        a(j + 1) = x
    end do
   end subroutine



