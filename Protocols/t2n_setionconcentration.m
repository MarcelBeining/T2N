function neuron = t2n_setionconcentration(neuron,str)


switch str
    case 'Krueppel'
        for t = 1:numel(neuron.mech)
            neuron.mech{t}.all.k_ion = struct('ek',-90);%-96.5);
            neuron.mech{t}.all.na_ion = struct('ena',60);%82.59);
            neuron.mech{t}.all.ca_ion = struct('cao0',2,'cai0',0.000048); % set hier also eca?
        end
    case 'Mongiat'
        for t = 1:numel(neuron.mech)
            neuron.mech{t}.all.k_ion = struct('ek',-93);%-100);%-103); !!-85phys!!!!!!   changed to -95 because there was lot of K outflux during current steps   % calculated from concentrations in Mongiat paper (SH07 very similar)
            neuron.mech{t}.all.na_ion = struct('ena',87.76);   %87 !63phys!!!!!!!!    % calculated from concentrations in Mongiat paper (SH07 very similar)
            neuron.mech{t}.all.ca_ion = struct('cao0',2,'cai0',0.000048);   %!!138phys!!!!!!    % calculated from concentrations in Mongiat paper (SH07 very similar)
        end
    case 'Mongiat2'  % if nakbuffer is available too...
        for t = 1:numel(neuron.mech)
            neuron.mech{t}.all.k_ion = struct('ko0',4,'ki0',140);%'ek',-93);%-100);%-103); !!-85phys!!!!!!   changed to -95 because there was lot of K outflux during current steps   % calculated from concentrations in Mongiat paper (SH07 very similar)
            neuron.mech{t}.all.na_ion = struct('nao0',156,'nai0',5);%87.76);   %87 !63phys!!!!!!!!    % calculated from concentrations in Mongiat paper (SH07 very similar)
            neuron.mech{t}.all.ca_ion = struct('cao0',2,'cai0',0.000048);   %!!138phys!!!!!!    % calculated from concentrations in Mongiat paper (SH07 very similar)
        end
    case {'Riazanski','SH07'}
        for t = 1:numel(neuron.mech)
            neuron.mech{t}.all.k_ion = struct('ek',-103.2);%-103); 
            neuron.mech{t}.all.na_ion = struct('ena',82.59);   
            neuron.mech{t}.all.ca_ion = struct('cao0',2,'cai0',0.000048);   
        end
    case 'SH08y' % action potential iniation and propagation in hc mossy fiber axons
        for t = 1:numel(neuron.mech)
            neuron.mech{t}.all.k_ion = struct('ek',-103.2);%-103); 
            neuron.mech{t}.all.na_ion = struct('ena',81.07);   
            neuron.mech{t}.all.ca_ion = struct('cao0',2,'cai0',0.000048);   
        end
    case 'SH08a'
        for t = 1:numel(neuron.mech)
            neuron.mech{t}.all.k_ion = struct('ek',-103.2);%-103); 
            neuron.mech{t}.all.na_ion = struct('ena',78.78);   
            neuron.mech{t}.all.ca_ion = struct('cao0',0.5,'cai0',0.000048);   
        end
end