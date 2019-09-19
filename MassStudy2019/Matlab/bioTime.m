% This code has timestamp information for all the BIOPAC data
% May 02, 2019

function dateVector = bioTime(fileName)
% Input: fileName: Name of biopac file
% Output: dateTime corresponding to input file in [Y M D H MI S] format

switch(fileName)
    case 'case1_biopac'
        dateVector = [2019,04,14,16,09,01];
        
    case 'case_add1_biopac'
        dateVector = [2019 04 20 12 20 40];
        
    case 'case2_biopac'
        dateVector = [2019 04 21 15 13 07];
        
    case 'case3_biopac'
        dateVector = [2019 04 22 20 16 55];
        
    case 'case4_biopac_1'
        dateVector = [2019 04 27 10 26 24];
        
    case 'case4_biopac_2'
        dateVector = [2019 04 27 12 23 34];
        
    case 'case_add2_13_06_biopac'
        dateVector = [2019 06 13 20 19 53];
        
    case 'case_add2_20_06_biopac'
        dateVector = [2019 06 20 14 47 36];
        
    case 'case_add2_29_06_biopac'
        dateVector = [2019 06 29 14 36 30];
        
    case 'case5_biopac'
        dateVector = [2019 06 18 14 56 57];
        
    case 'case_add2_BT_07_01_biopac'
        dateVector = [2019 07 01 15 16 56];
        
    case 'case_add2_BT_07_04_biopac'
        dateVector = [2019 07 04 13 29 16];
        
    case 'case_add3_BT_07_03_biopac'
        dateVector = [2019 07 03 16 06 42];
        
    otherwise
        fprintf('Enter valid BIOPAC filename\n');
end

end