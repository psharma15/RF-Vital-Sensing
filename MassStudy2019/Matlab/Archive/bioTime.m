% This code has timestamp information for all the BIOPAC data
% May 02, 2019

function dateVector = bioTime(fileName)
% Input: fileName: Name of biopac file
% Output: dateTime corresponding to input file in [Y M D H MI S] format

switch(fileName)
    case 'Case1_14April_biopac'
        dateVector = [2019,04,14,16,09,01];
        
    otherwise
        fprintf('Enter valid BIOPAC filename\n');
end

end