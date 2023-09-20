
function [cycle_scat,cycle_pat,cycle_neg,pat_area,Initial_pos,Final_pos,Initial_AI,Final_AI] = OPaL_Single_Run(size,nuc_wid,time,iterations,threshold,induction,ca,aa)

presence = 0;
presence_threashold = 4;
scat = 0.5;
pat = 0.4;


 
%presence coutner for angle
pres_count = 0;
pres_count_initial = 0;
pres_count_final = 0;
pres_incr = 0;
jj = 0;
buff_angle = 0;
one_and_one = 0; %indicator of 1 at the beggining and end
new_jj = 1;
new_jjstart = 1;
 


 
%calculating the convergence of 
converge(iterations,time)=0;
converge_count = 1;

cycle_pat = 0;
cycle_neg = 0;
cycle_scat = 0;

row = 0;
col = 0;
        
Pos_matrix(time,size) = 0;
AI_matrix(time,size) = 0;
Temp_AI(size,size) = 0;
Final_pos(iterations,size) = 0;
Initial_pos(iterations,size) = 0;
Initial_AI(iterations,size) = 0;
Final_AI(iterations,size) = 0;
pat_area(iterations,10)=0;
excel(iterations*2,size)=0;
 

       
 

%position value, activiation matrix
for i = 1:size
    if i == 1
        for j = 1:size
            if j <= ((size/2)+1)
                xx(i,j)= j;
            else
                xx(i,j)= xx(1,j-1)-1;
            end
            A_matrix(i,j) = aa*exp(-ca*xx(i,j)*nuc_wid);
            I_matrix(i,j) = 0;
        end
    else
        for j = 1:size
            if j == 1
                xx(i,j) = xx(i-1,size);
            else
                xx(i,j) = xx(i-1,j-1);               
            end
            A_matrix(i,j) = aa*exp(-ca*xx(i,j)*nuc_wid);
            I_matrix(i,j) = 0;
        end
    end
end


% Position matrix
for ite = 1:iterations
    for i = 1:time
        if i == 1
            for j = 1:size
                Pos_matrix(i,j) = randi([0,1]);                            
                if rand <= induction
                    Pos_matrix(i,j) = 1;
                else
                    Pos_matrix(i,j) = 0;
                end
            end
            %Pos_matrix(i,18:29) = 1;
            for j = 1:size
                if Pos_matrix(i,j) > 0
                    Temp_AI = A_matrix.*transpose(Pos_matrix(i,:))-I_matrix.*transpose(Pos_matrix(i,:));
                end
            end
            for j = 1:size
                AI_matrix(i,:) = sum(Temp_AI);
            end
        else
            for j = 1:size
                if AI_matrix(i-1,j) >= threshold
                    Pos_matrix(i,j) = 1;
                else
                    Pos_matrix(i,j) = 0;
                end
            end

            Temp_AI = A_matrix.*transpose(Pos_matrix(i,:))-I_matrix.*transpose(Pos_matrix(i,:));
            for j = 1:size
                AI_matrix(i,:) = sum(Temp_AI);
            end
        end
    end
    for i = 1:time
       converge_value = sum(abs(Pos_matrix(time,:)-Pos_matrix(i,:)));
       converge(converge_count,i)= converge_value;
    end
    converge_count = converge_count + 1;
    Initial_pos(ite,:)= Pos_matrix(1,:);
    Final_pos(ite,:)= Pos_matrix(time,:);
    Initial_AI(ite,:) = AI_matrix(1,:);
    Final_AI(ite,:) = AI_matrix(time,:);
    excel(ite*2-1,:) = Pos_matrix(1,:);
    excel(ite*2,:) = Pos_matrix(time,:);
end

sum_positions = sum(Final_pos,2);
converge = 1-(converge/size);

%patterning evaluation
for i = 1:iterations
    presence = 0;
    pres_incr = 0;
    pres_count_initial = 0;
    pres_count_final = 0;
    new_jj = size;
    new_jjstart = 1;
    if Final_pos(i,1) == 1 && Final_pos(i,size) == 1
        pres_incr = 0;
        one_and_one = 1;
        for jj = 1:(size-1)
            if pres_incr == 0
                if abs(Final_pos(i,jj)-Final_pos(i,jj+1))>0
                    pres_incr = pres_incr + 1;
                    pres_count_initial = jj;
                    new_jjstart = jj+2;
                end
            end
        end
        jj = size;
        for jj = size:-1:1
            if pres_incr == 1
                if abs(Final_pos(i,jj)-Final_pos(i,jj-1))>0
                    pres_incr = pres_incr + 1;
                    pres_count_final = size-jj;
                    buff_angle = (360/size)*(((pres_count_final+pres_count_initial)/2)+jj);
                    new_jj = jj-1;
                    if buff_angle > 360
                        pat_area(i,4) = buff_angle-360;
                    else
                        pat_area(i,4) = buff_angle;
                    end
                end
            end
        end
        presence = 2;
        pres_incr = 0;
    end

    if Final_pos(i,1) == 0 && Final_pos(i,size) == 1
        pres_incr = 0;
        pres_count_initial = size;
        pres_incr = pres_incr + 1;
        for jj = size:-1:1
            if pres_incr == 1
                if abs(Final_pos(i,jj)-Final_pos(i,jj-1))>0
                    pres_incr = pres_incr + 1;
                    pres_count_final = jj;
                    buff_angle = (360/size)*(((pres_count_final+pres_count_initial)/2));
                    new_jj = jj-1;
                    if buff_angle > 360
                        pat_area(i,4) = buff_angle-360;
                    else
                        pat_area(i,4) = buff_angle;
                    end
                end
            end
        end
        presence = 2;
        pres_incr = 0;
    end

    if Final_pos(i,1) == 1 && Final_pos(i,size) == 0
        pres_incr = 0;
        pres_count_final = 0;
        pres_incr = pres_incr + 1;
        for jj = 1:(size-1)
            if pres_incr == 1
                if abs(Final_pos(i,jj)-Final_pos(i,jj+1))>0
                    pres_incr = pres_incr + 1;
                    pres_count_initial = size-jj;
                    buff_angle = (360/size)*(((pres_count_final+pres_count_initial)/2)+jj);
                    new_jjstart = jj+2;
                    if buff_angle > 360
                        pat_area(i,4) = buff_angle-360;
                    else
                        pat_area(i,4) = buff_angle;
                    end
                end
            end
        end
        presence = 2;
        pres_incr = 0;
    end

    for j = new_jjstart:1:new_jj-1
        if abs(Final_pos(i,j)-Final_pos(i,j+1))>0
            presence = presence + 1;
        end
        if abs(Final_pos(i,j)-Final_pos(i,j+1))>0
            pres_incr = pres_incr + 1;
            if pres_incr == 1
                pres_count_initial = j+1;
            end
            if pres_incr == 2
                pres_count_final = j+1;
                if pres_count_initial == 1
                    pat_area(i,3+(presence/2))= (360/size)*((pres_count_final - pres_count_initial+1)/2+j);

                else
                    pat_area(i,3+(presence/2))= (360/size)*((pres_count_final - pres_count_initial+1)/2+(pres_count_initial-1));
                end
                pres_incr = 0;
            end                             
        end
    end

    if sum_positions(i)/size == 0
         pat_area(i,1) = -1;
    else
        if (sum_positions(i)/size > scat) || ((sum_positions(i)/size > pat) && (presence/2 > presence_threashold))
            pat_area(i,1) = 0;
        else
            pat_area(i,1) = 1;
        end
    end
    
    pat_area(i,2) = sum_positions(i)/size;
    pat_area(i,3) = presence/2;
    pres_incr = 0;
    one_and_one = 0;
    
end
cycle_pat = cycle_pat + sum(pat_area(:,1)==1)/iterations;
cycle_neg = cycle_neg + sum(pat_area(:,1)==-1)/iterations;
cycle_scat = cycle_scat + sum(pat_area(:,1)==0)/iterations;