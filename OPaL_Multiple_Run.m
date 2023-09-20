
function [result_matrix_scat,result_matrix,result_matrix_neg] = OPaL_Multiple_Run(size,nuc_wid,time,iterations,induction,aa,threshold_ini,threshold_end,threshold_inc,ca_ini,ca_end,ca_inc)
    presence = 0;
    presence_threashold = 4;
    scat = 0.5;
    pat = 0.4;
    reps = 1;


    % For high throuput
    cycle_pat = 0;
    cycle_neg = 0;
    cycle_scat = 0;

    row = 0;
    col = 0;


    result_matrix(fix((threshold_end-threshold_ini)/threshold_inc),fix((ca_end-ca_ini)/ca_inc)) = 0;
    result_matrix_neg(fix((threshold_end-threshold_ini)/threshold_inc),fix((ca_end-ca_ini)/ca_inc)) = 0;
    result_matrix_scat(fix((threshold_end-threshold_ini)/threshold_inc),fix((ca_end-ca_ini)/ca_inc)) = 0;


    for threshold = threshold_ini:threshold_inc:threshold_end
        row = row + 1;
        col = 0;

        for ca = ca_ini:ca_inc:ca_end
                   col = col+1;

            xx = 0;
            A_matrix = 0;
            I_matrix = 0;
            for k = 1:reps
                Pos_matrix(time,size) = 0;
                AI_matrix(time,size) = 0;
                Temp_AI(size,size) = 0;
                Final_pos(iterations,size) = 0;
                Initial_pos(iterations,size) = 0;
                sum_positions(iterations, size) = 0;
                pat_area(iterations,3)=0;
                excel(iterations*2,size)=0;



                if k == 1
                % position value, activiation and inhibition matrix
                    for i = 1:size
                        if i == 1
                            for j = 1:size
                                if j <= ((size/2)+1)
                                    xx(i,j)= j;
                                else
                                    xx(i,j)= xx(1,j-1)-1;
                                end
                                A_matrix(i,j) = aa*exp(-ca*xx(i,j)*nuc_wid);
                                %I_matrix(i,j) = ai*exp(-ci*xx(i,j)*nuc_wid);
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
                                %I_matrix(i,j) = ai*exp(-ci*xx(i,j)*nuc_wid);
                                I_matrix(i,j) = 0;
                            end
                        end
                    end
                end

                % Position matrix
                for ite = 1:iterations
                    for i = 1:time
                        if i == 1
                            for j = 1:size
                               % Pos_matrix(i,j) = randi([0,1]);

                                if rand <= induction
                                    Pos_matrix(i,j) = 1;
                                else
                                    Pos_matrix(i,j) = 0;
                                end

                            end
                    %         for j = 1:size
                    %             if Pos_matrix(i,j) > 0
                                    Temp_AI = A_matrix.*transpose(Pos_matrix(i,:))-I_matrix.*transpose(Pos_matrix(i,:));
                    %             end
                    %         end
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
                    Initial_pos(ite,:)= Pos_matrix(1,:);
                    Final_pos(ite,:)= Pos_matrix(time,:);  
                    excel(ite*2-1,:) = Pos_matrix(1,:);
                    excel(ite*2,:) = Pos_matrix(time,:);
                end

                sum_positions = sum(Final_pos,2);

                %patterning evaluation

                for i = 1:iterations
                    presence = 0;
                    for j = 1:size
                        if j == 1
                            if abs(Final_pos(i,size)-Final_pos(i,j))>0
                                presence = presence + 1;
                            end
                        else
                            if abs(Final_pos(i,j)-Final_pos(i,j-1))>0
                                presence = presence + 1;
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
                end
                cycle_pat = cycle_pat + sum(pat_area(:,1)==1)/iterations;
                cycle_neg = cycle_neg + sum(pat_area(:,1)==-1)/iterations;
                cycle_scat = cycle_scat + sum(pat_area(:,1)==0)/iterations;
                clear Pos_matrix;
                clear AI_matrix;
                clear Temp_AI;
                clear Final_pos;
                clear Initial_pos;
                clear sum_positions;
                clear pat_area;
                clear excel;


            end
            result_matrix(row,col) = cycle_pat/reps;
            result_matrix_neg(row,col) = cycle_neg/reps;
            result_matrix_scat(row,col) = cycle_scat/reps;
            cycle_pat = 0;
            cycle_neg = 0;
            cycle_scat = 0;
            clear xx;
            clear A_matrix;
            clear I_matrix;
        end
    end



