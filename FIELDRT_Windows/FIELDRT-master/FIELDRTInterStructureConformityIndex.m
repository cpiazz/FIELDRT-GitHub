function output = FIELDRTInterStructureConformityIndex(filelist,pathname,struct_name,refStudyName,type)
%
% Calculate conformity index for a structure among a set of CERR containers.
% Filelist must be a cell array (as from uigetfile).
%
% Types of conformity index:
%
% type 1 
% ----------
% RTOG conformity index CI = Vri / TV
% 
% Vri = volume covered by reference isodose
% TV = target volume
%
% type 2
% ----------
% van't Riet et al. conformity number CN = (TVri / TV) * (TVri / Vri)
% 
% Vri = volume covered by reference isodose (structure volume)
% TVri = target volume covered by reference isodose (intersection volume)
% TV = target volume (gold standard structure volume)
%
% type 3
% ----------
% Jaccard similarity coefficient CI = (Vri /\ TV) / (Vri \/ TV)
% 
% Vri = volume covered by reference isodose (structure volume)
% TV = target volume (gold standard structure volume)
% /\ = intersection operator
% \/ = union operator
%
% type 4
% ----------
% Dice coefficient CI = 2 * (Vri /\ TV) / (Vri + TV)
% 
% Vri = volume covered by reference isodose (structure volume)
% TV = target volume (gold standard structure volume)
% /\ = intersection operator
%
% type 5
% ----------
% Geographical Miss Index GMI = [TV - (Vri /\ TV)] / TV
%
% Vri = volume covered by reference isodose (structure volume)
% TV = target volume (gold standard structure volume)
% /\ = intersection operator
% \/ = union operator
%
% type 6
% ----------
% Discordance Index DI = 1 - (Vri /\ TV)/Vri 
%
% Vri = volume covered by reference isodose (structure volume)
% TV = target volume (gold standard structure volume)
% /\ = intersection operator
% \/ = union operator

% Concetta Piazzese: last update MAR 2018


% nfiles=length(filelist);
% output=cell(nfiles,2);

[refVol,structRefNum,planCref]=getRefVol(filelist,pathname,refStudyName,struct_name);

% for k=1:nfiles
%     filename=fullfile(pathname,filelist{k});
%     disp(['Processing file: ',filename]);
%     load(filename);
%     planCeval=planC;
%     clear planC
    indexS=planCeval{end};
    % Get file name
    output{k,1}=filelist{k};
    % Get structure number
    structEvalNum=getStructNum(struct_name,planCeval,indexS);
    % Get structure volume
    vol=getfield(structureStats(structEvalNum,planCeval), 'vol');
    % Get intersection/union volumes
    [unionVolume,intersectVolume]=calcInterUnionIntersectionVol(structEvalNum,structRefNum,planCeval,planCref);
%     [unionVolume,intersectVolume]=calcInterUnionIntersectionArea(structEvalNum,structRefNum,planCeval,planCref);
    
    % Calculate and save conformity index
%     switch type
%         case 1 % RTOG
%             output{k,2}=vol/refVol;
%         case 2 % van't Riet
%             output{k,2}=(intersectVolume/refVol)*(intersectVolume/vol);            
%         case 3 % Jaccard
%             output{k,2}=intersectVolume/unionVolume;
%         case 4 % Dice
%             output{k,2}=2*intersectVolume/(vol+refVol);
%     end
    switch type
        case 'RTOG' % 
            output{k,2}=vol/refVol;
        case 'RIET' % 
            output{k,2}=(intersectVolume/refVol)*(intersectVolume/vol);            
        case 'JACCARD' % 
            output{k,2}=intersectVolume/unionVolume;
        case 'DICE' % 
            output{k,2}=2*intersectVolume/(vol+refVol);
        case 'GMI' % 
            output{k,2}=(refVol-intersectVolume)/refVol;
        case 'DI' % 
            output{k,2}=1-(intersectVolume/vol);
    end
% end
toc


function [refVol,refStructNum,planC]=getRefVol(filelist,pathname,refStudyName,struct_name)

% Find the number of the dataset in the list which corresponds to the reference dataset 
index = strfind(filelist, refStudyName);

% Retrieve reference volume and reference structure number in reference
% dataset. Multiple entries will be ignored. No error check is made.
for k=1:length(index)
    if index{k}==1
        filename=fullfile(pathname,filelist{k});
        load(filename);
        indexS = planC{end};
        numStruct=length(planC{indexS.structures});
        refStructNum=getStructNum(struct_name,planC,indexS);
        refVol=getfield(structureStats(refStructNum,planC), 'vol');
    end
end        
