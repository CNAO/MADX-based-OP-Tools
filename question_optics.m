% {}~
% - include Matlab library
pathToLibrary=".\externals\MatLabTools";
addpath(genpath(pathToLibrary));
% OPpath="K:";
OPpath="S:\Accelerating-System\Accelerator-data";

% magnetName="V1_044A_CEB"; %[ "V1_044A_CEB" ; "V2_013A_CEB" ]; % [ "T1_011A_CEB" "T2_015A_CEB" ],[ "U1_023A_CEB" "U2_013A_CEB" ],[ "V1_044A_CEB" "V2_013A_CEB" ],[ "Z1_011A_CEB" "Z2_015A_CEB" ]
% properties=[ ... % 1 row per magnet
%     "HKICK" "VKICK" ; ... % first corrector
%     % "HKICK" "VKICK" ; ... % second corrector
%     ];
% observations=[ ... % observations in MADX output files, will be used in cartesian product with magnet properties
%     "X:V2_031A_MOB" ...
%     "Y:V2_031A_MOB" ...
%     "X:V2_029B_NZL" ...
%     "Y:V2_029B_NZL" ...
%     "X:V2_019B_SFH" ...
%     "Y:V2_019B_SFH" ...
%     ];
% for iMagnet=1:size(magnetNames,1)
%     whatToShow(:,:,iMagnet)=[ ...
%         "HKICK:X:V2_031A_MOB" "HKICK:X:V2_029B_NZL" "HKICK:X:V2_019B_SFH" ; ...
%         "VKICK:Y:V2_031A_MOB" "VKICK:Y:V2_029B_NZL" "VKICK:Y:V2_019B_SFH" ; ...
% %         "HKICK:X@V2_031A_MOB" ; ...
%     ];
%     fitOrders(:,:,iMagnet)=[ ...
%         1 1 1 ; ...
%         1 1 1 ; ...
%     ];
% end
% % some checks
% if (size(magnetNames,1)~=size(properties,1) )
%     error("Not all magnets are assigned a property to vary!");
% end
% if (size(magnetNames,1)~=size(whatToShow,3) )
%     error("Not all magnets are assigned a visualisation scheme!");
% end
% % check that all visualisation schemes are found in the data parsed -- to come!

beamPart="CARBON";   % select beam particle: PROTON, CARBON
outputPath="OUTPUT"; % for .xlsx files
MADXpath=".\externals\optics\HEBT";        % path to MADX files
combo="REGULAR";     % REGULAR, COUPLED -- see header of getStandards function
mode="RM";           % RM: response matrix
% magnetName="T1_011A_CEB";
% magnetName="T2_015A_CEB";
magnetName="U1_023A_CEB";
% magnetName="U2_013A_CEB";
% magnetName="V1_044A_CEB";
% magnetName="V2_013A_CEB";
% magnetName="Z1_011A_CEB";
% magnetName="Z2_015A_CEB";
Is=[ -10 10 10 ];   % Imin, Imax, Istep [A]
NeNs=[ 1 -1 1 ];      % energy values: min, max, step -- if max<min, go up to last energy value
% properties=[ "HKICK" "VKICK" ];
% pippo=[ ... % observations in MADX output files
%     "X:V2_031A_MOB" ...
%     "Y:V2_031A_MOB" ...
%     "X:V2_029B_NZL" ...
%     "Y:V2_029B_NZL" ...
%     "X:V2_019B_SFH" ...
%     "Y:V2_019B_SFH" ...
%     ];
% whatToShow=[ ...
%     "HKICK:X:V2_031A_MOB" "HKICK:X:V2_029B_NZL" "HKICK:X:V2_019B_SFH" ; ...
%     "VKICK:Y:V2_031A_MOB" "VKICK:Y:V2_029B_NZL" "VKICK:Y:V2_019B_SFH" ; ...
% ];
% fitOrders=[ ...
%     1 1 1 ; ...
%     1 1 1 ; ...
% ];

% pre-processing
beamPart=upper(beamPart);
magnetName=upper(magnetName);
% standard plots, properties and observations
[whatToShow,fitOrders]=getStandards(magnetName,mode,combo);
[properties,observations]=WTS_2_PROP_OBS(whatToShow);
% standard filenames
[xlsFileName,MADXfileNames,RMfileNames]=fileNames(magnetName,properties);

writeMADXfiles(magnetName,properties,observations,beamPart,MADXpath,MADXfileNames,RMfileNames,Is,NeNs); % generate MADX files
runMADX(MADXfileNames(2),MADXpath); % run MADX
[RMtables,RMtableHeaders]=readMADXData(MADXpath,RMfileNames);
overview3D(RMtables,RMtableHeaders,magnetName,whatToShow,properties);
[FitParams,headerFitParams]=performFit(RMtables,RMtableHeaders,magnetName,whatToShow,properties,fitOrders);
writeXLStable(FitParams,headerFitParams,xlsFileName,outputPath,whatToShow,1);
% [FitParams,headerFitParams]=readXLStable(xlsFileName,outputPath,whatToShow);
% showFitParams(FitParams,magnetName,whatToShow);
showFitParams(FitParams,magnetName,whatToShow,OPpath,beamPart);
showFitParams(FitParams,magnetName,whatToShow,OPpath,beamPart,1); % show ratio wrt OP RMs

%% MADX

function writeMADXfiles(magnetName,properties,observations,beamPart,MADXpath,MADXfileNames,RMfileNames,Is,NeNs)
    templateMasterFile=MADXfileNames(1);
    actualMasterFile=MADXfileNames(2);
    templateRMFile=MADXfileNames(3);
    actualRMFile=MADXfileNames(4);
    origFolder=cd(MADXpath);
    
    fprintf('generating MADX master file %s out of %s ...\n',actualMasterFile,templateMasterFile);
    % read in template file
    fid = fopen(templateMasterFile, 'r');
    C=textscan(fid,'%s','delimiter','\n'); % this is a cell array
    fclose(fid);
    % modify concerned lines
    templateLen=-1;
    for k=1:numel(C{1,1})
        if ( startsWith(strip(string(C{1,1}(k))),"is_carbon=") )
            switch beamPart
                case "PROTON"
                    is_carbon=0;
                case "CARBON"
                    is_carbon=1;
                otherwise
                    error("Unable to identify beam particle %s!",beamPart);
            end
            C{1,1}(k)=cellstr(sprintf("is_carbon=%.0f;",is_carbon));
        elseif ( startsWith(strip(string(C{1,1}(k))),"iLine=") )
            switch extractBetween(magnetName,1,1)
                case "T"
                    iLine=1;
                case "U"
                    iLine=2;
                case "V"
                    iLine=3;
                case "Z"
                    iLine=4;
                otherwise
                    error("Unable to identify line from magnet name %s!",magnetName);
            end
            C{1,1}(k)=cellstr(sprintf("iLine=%.0f;",iLine));
        elseif ( startsWith(strip(string(C{1,1}(k))),'! end of template') )
            templateLen=k;
            break
        end
    end
    % write it out
    fout = fopen(actualMasterFile, 'w');
    for k=1:templateLen 
        fprintf(fid,'%s\r\n',C{1,1}{k,1}); 
    end
    fprintf(fid,sprintf('call, file="%s";',actualRMFile));
    fclose(fout);

    fprintf('generating MADX RM file %s out of %s ...\n',actualRMFile,templateRMFile);
    % read in template file
    fid = fopen(templateRMFile, 'r');
    C=textscan(fid,'%s','delimiter','\n'); % this is a cell array
    fclose(fid);
    % modify concerned lines
    templateLen=-1;
    for k=1:numel(C{1,1})
        if ( startsWith(strip(string(C{1,1}(k))),"Imin=") )
            C{1,1}(k)=cellstr(sprintf("Imin=%.3f; ! [A]",Is(1)));
        elseif ( startsWith(strip(string(C{1,1}(k))),"Imax=") )
            C{1,1}(k)=cellstr(sprintf("Imax=%.3f; ! [A]",Is(2)));
        elseif ( startsWith(strip(string(C{1,1}(k))),"Istep=") )
            C{1,1}(k)=cellstr(sprintf("Istep=%.3f; ! [A]",Is(3)));
        elseif ( startsWith(strip(string(C{1,1}(k))),"NEnLevelstart=") )
            C{1,1}(k)=cellstr(sprintf("NEnLevelstart=%.0f;",NeNs(1)));
        elseif ( startsWith(strip(string(C{1,1}(k))),"NEnLevelstop=") )
            if ( NeNs(2)<NeNs(1) )
                C{1,1}(k)=cellstr("NEnLevelstop=lunghezzaLGEN;");
            else
                C{1,1}(k)=cellstr(sprintf("NEnLevelstop=%.0f;",NeNs(2)));
            end
        elseif ( startsWith(strip(string(C{1,1}(k))),"NEnLevelstep=") )
            C{1,1}(k)=cellstr(sprintf("NEnLevelstep=%.0f;",NeNs(3)));
        elseif ( startsWith(strip(string(C{1,1}(k))),'! end of template') )
            templateLen=k;
            break
        end
    end
    % write it out
    fid = fopen(actualRMFile, 'w');
    for k=1:templateLen 
        fprintf(fid,'%s\r\n',C{1,1}{k,1}); 
    end
    % write macros
    % - write header macro
    fprintf(fid,'! header of output file\r\n');
    fprintf(fid,'writeRMHeader(fileName): macro{\r\n');
    fprintf(fid,'   assign, echo=fileName;\r\n');
    [iCols,labels,units]=fixedColsMADXfile();
    fixedCols=strings(length(iCols),1);
    for ii=1:length(iCols)
        [~,fixedCols(ii)]=getLabels(labels(ii),iCols,labels,units);
    end
    fixString=sprintf(" %s,",fixedCols);
    fixString=extractBetween(fixString,1,strlength(fixString)-1);
    % add units
    tmpObs=observations;
    for ii=1:length(tmpObs)
        tmpSplit=split(tmpObs(ii),":");
        tmpObs(ii)=sprintf("%s[%s]",tmpObs(ii),MADXunits(tmpSplit(1)));
    end
    obsString=sprintf(" %s,",tmpObs);
    obsString=extractBetween(obsString,1,strlength(obsString)-1);
    fprintf(fid,'   PRINT, TEXT="# %s, %s";\r\n',fixString,obsString);
    fprintf(fid,'   assign, echo=terminal;\r\n');
    fprintf(fid,'};\r\n');
    % - write computing macro
    fprintf(fid,'! actual macro\r\n');
    fprintf(fid,'writeRMObservations(fileName): macro{\r\n');
    fprintf(fid,'   assign, echo=fileName;\r\n');
    obsColsFormat=strings(length(observations),1);
    obsColsFormat(:)="% 19.12E ";
    obsColsFormat=sprintf("%s,",obsColsFormat);
    obsColsFormat=extractBetween(obsColsFormat,1,strlength(obsColsFormat)-1);
    fprintf(fid,'   PRINTF, TEXT="%% 19.12E,%% 6.1f,%% 12.5E,%% 19.12E,%s",\r\n',obsColsFormat);
    fprintf(fid,'           VALUE=Brho, BP, currI, kWrite,\r\n');
    for ii=1:length(observations)
        allSplit=split(observations,":");
        if ( ii==length(observations) )
            lastChar=';'; % last char closes MADX command!
        else
            lastChar=',';
        end
        obString=sprintf('           table(twiss,%s,%s)%s\r\n',allSplit(1,ii,2),allSplit(1,ii,1),lastChar);
        fprintf(fid,obString);
    end
    fprintf(fid,'   assign, echo=terminal;\r\n');
    fprintf(fid,'};\r\n');
    for iProperty=1:length(properties)
        fprintf(fid,'exec, compute_RM_allCurrents(%s,%s,I2K_HEBT_CEBH,%s);\r\n',magnetName,properties(iProperty),RMfileNames(iProperty));
    end
    fclose(fout);
    
    cd(origFolder);
end

% run MADX
function runMADX(actualMasterFile,MADXpath)
    origFolder=cd(MADXpath);
    MADexe="T:\ARicerca\MADX\5.06.01\madx-win64-gnu.exe";
    % MADexe="K:\Area dati MD\00MatriciDiRisposta\MADX-based-OP-Tools\madx5.06.01\madx-win64-gnu.exe";
    fprintf('running %s ...\n',MADexe);
    % mycmd=sprintf('"%s" < %s > MADX.log',MADexe,actualMasterFile);
    mycmd=sprintf('"%s" < %s',MADexe,actualMasterFile);
    [status,cmdout]=system(mycmd);
    cd(origFolder);
    if ( status ~= 0 )
        error("error in running MADX!");
    end
end

% read data files
function [RMtables,RMtableHeaders]=readMADXData(MADXpath,RMfileNames)
    origFolder=cd(MADXpath);
    [iCols,~,~]=fixedColsMADXfile();
    nFixedColumns=length(iCols);
    nFiles=length(RMfileNames);

    % parsing MADX data
    nMaxCols=nFixedColumns;
    nMaxRows=1;
    RMtables=zeros(nMaxRows,nMaxCols,nFiles);   % structure of RMtables
    RMtableHeaders=strings(1,nMaxCols,nFiles);  % structure of RMtableHeaders
    for iFile=1:nFiles
        % MADX file name and path:
        fprintf('parsing file %s ...\n',RMfileNames(iFile));
        RMtable=readmatrix(RMfileNames(iFile),'Delimiter',',','NumHeaderLines',1,'FileType','text');
        % increase storage, if needed
        if ( size(RMtable,1)>size(RMtables,1) | size(RMtable,2)>size(RMtables,2) )
            nMaxRows=max(size(RMtable,1),size(RMtables,1));
            nMaxCols=max(size(RMtable,2),size(RMtables,2));
            % data table
            temp=RMtables;
            RMtables=zeros(nMaxRows,nMaxCols,nFiles);
            RMtables(1:size(temp,1),1:size(temp,2),1:nFiles)=temp;
            clear temp;
            % header
            temp=RMtableHeaders;
            RMtableHeaders=strings(1,nMaxCols,nFiles);
            RMtableHeaders(1,1:size(temp,2),1:nFiles)=temp;
            clear temp;
        end
        RMtables(1:size(RMtable,1),1:size(RMtable,2),iFile)=RMtable;
        clear RMtable;
        % get header
        fid = fopen(RMfileNames(iFile), 'r');
        header = fgets(fid);
        fclose(fid);
        temp=strip(split(string(header(2:end)),","))';
        RMtableHeaders(1,1:size(temp,2),iFile)=temp;
    end
    cd(origFolder);
end

% fixed format of data files
function [iCols,labels,units]=fixedColsMADXfile()
    iCols=1:4;
    labels=[ "Brho" "BP" "I" "K" ];
    units=[ "Tm" "mm" "A" "rad" ];
end

%% specific functions
function overview3D(RMtables,RMtableHeaders,magnetName,whatToShow,properties)
%   make 3D plots giving an overview of the dependence
    fprintf("3D data overview...\n");
    [iCols,labels,units]=fixedColsMADXfile(); 
    xVarName="I";  [iColX,labelX]=getLabels(xVarName,iCols,labels,units);
    parName="BP";  [iColPar,labelPar]=getLabels(parName,iCols,labels,units);
    myTitle=LabelMe(magnetName);
    currTitle=sprintf("%s - 3D overview",myTitle);
    ff=figure('Name',currTitle,'NumberTitle','off');
    % get unique set of BP values
    sets=unique(RMtables(:,iColPar,:));
    nSets=length(sets);
    nRows=size(whatToShow,1);
    nCols=size(whatToShow,2);
    cm=colormap(parula(nSets));
    for iRow=1:nRows
        for iCol=1:nCols
            iPlot=(iRow-1)*nCols+iCol;
            subplot(nRows,nCols,iPlot);
            tempSplit=split(whatToShow(iRow,iCol),":");
            kick=tempSplit(1);
            obsVar=tempSplit(2);
            obsLoc=tempSplit(3);
            obs=sprintf("%s:%s",obsVar,obsLoc);
            propIndex=find(properties==kick);
            obsIndex=contains(RMtableHeaders(:,:,propIndex),obs);
            [multFact,unit]=FactorMe(obsVar);
            for iSet=1:nSets
                indices=(RMtables(:,iColPar,propIndex)==sets(iSet));
                plot3(RMtables(indices,iColX,propIndex),...
                      RMtables(indices,iColPar,propIndex),...
                      RMtables(indices,obsIndex,propIndex)*multFact,...
                      'Color',cm(iSet,:));
                hold on;
            end
            title(sprintf("%s - %s",kick,LabelMe(obsLoc)));
            grid on;
            xlabel(labelX);
            ylabel(labelPar);
            zlabel(sprintf("%s [%s]",obsVar,unit));
        end
    end
    sgtitle(currTitle);
end

function [FitParams,headerFitParams]=performFit(RMtables,RMtableHeaders,magnetName,whatToShow,properties,fitOrders)
%   fit data
%   RM(I), parametric in BP
    fprintf("fitting data...\n");
    [iCols,labels,units]=fixedColsMADXfile(); 
    xVarName="I";  [iColX,labelX]=getLabels(xVarName,iCols,labels,units);
    parName="BP";  [iColPar,labelPar]=getLabels(parName,iCols,labels,units);
    [iColBrho,labelBrho]=getLabels("Brho",iCols,labels,units);
    myTitle=LabelMe(magnetName);
    currTitle=sprintf("%s - Fits",myTitle);
    ff=figure('Name',currTitle,'NumberTitle','off');
    % get unique set of BP values
    sets=unique(RMtables(:,iColPar,:));
    nSets=length(sets);
    nRows=size(whatToShow,1);
    nCols=size(whatToShow,2);
    nOrders=max(fitOrders,[],'all');
    % structure of FitParams
    % write BP,Brho,RM,p(:)'
    FitParams=zeros(nSets,nOrders+4,size(whatToShow,1),size(whatToShow,2));
    headerFitParams=strings(nOrders+4,1); % structure of header
    [iCols,labels,units]=fixedColsFitParams();
    for ii=1:length(iCols)
        [~,headerFitParams(ii)]=getLabels(labels(ii),iCols,labels,units);
    end
    iCol=length(iCols);
    for iOrder=nOrders:-1:2
        iCol=iCol+1;
        [~,headerFitParams(iCol)]=sprintf("p_%d [mm/A^%d]",iOrder,iOrder);
    end
    headerFitParams(end-1)="p_1 [mm/A]";
    headerFitParams(end)="p_0 [mm]";
    headerFitParams=headerFitParams';
    for iSet=1:nSets
        indices=(RMtables(:,iColPar,:)==sets(iSet));
        for iRow=1:nRows
            for iCol=1:nCols
                tempSplit=split(whatToShow(iRow,iCol),":");
                kick=tempSplit(1);
                obsVar=tempSplit(2);
                obsLoc=tempSplit(3);
                obs=sprintf("%s:%s",obsVar,obsLoc);
                propIndex=find(properties==kick);
                obsIndex=contains(RMtableHeaders(:,:,propIndex),obs);
                [multFact,unit]=FactorMe(obsVar);
                % indices=(RMtables(:,iColPar,propIndex)==sets(iSet));
                RMtempTab=RMtables(indices(:,:,propIndex),:,propIndex);
                % perform fit
                % xFit=RMtables(indices,iColX,propIndex);
                % yFit=RMtables(indices,obsIndex,propIndex)*multFact;
                xFit=RMtempTab(:,iColX);
                yFit=RMtempTab(:,obsIndex)*multFact;
                p = polyfit(xFit,yFit,fitOrders(iRow,iCol));
                xs=min(xFit):(max(xFit)-min(xFit))/200:max(xFit);
                ys = polyval(p,xs);
                % save fitting data
                % currBrho=unique(RMtables(indices,iColBrho,propIndex));
                currBrho=unique(RMtempTab(:,iColBrho));
                RM=p(end-1)*currBrho;
                FitParams(iSet,1:4+fitOrders(iRow,iCol),iRow,iCol)=[currBrho sets(iSet) RM p(:)' ];
                % plot
                iPlot=(iRow-1)*nCols+iCol;
                subplot(nRows,nCols,iPlot);
                plot(xFit,yFit,'k.',xs,ys,'r-');
                legend('data','fit','Location','best');
                grid on;
                xlabel(labelX);
                ylabel(sprintf("%s [%s]",obsVar,unit));
                title(sprintf("%s - %s",kick,LabelMe(obsLoc)));
            end
        end
        drawnow;
        sgtitle(sprintf("%s - %s=%.1f",currTitle,labelPar,sets(iSet)));
        % pause(0.0);
    end
end

function [iCols,labels,units]=fixedColsFitParams()
    iCols=1:3;
    labels=[ "Brho" "BP" "RM" ];
    units=[ "Tm" "mm" "mm/A Tm"  ];
end

function showFitParams(FitParams,magnetName,whatToShow,OPpath,beamPart,iRatio)
    fprintf("showing dependence of RMs...\n");
    [iCols,labels,units]=fixedColsFitParams();
    xVarName="BP";     [iColX,labelX]=getLabels(xVarName,iCols,labels,units); % used for plotting
    yVarName="RM";     [iColY,labelY]=getLabels(yVarName,iCols,labels,units); % used for plotting
    xVarNameRM="Brho"; [iColXrm,labelXrm]=getLabels(xVarNameRM,iCols,labels,units);
    myTitle=LabelMe(magnetName);
    currTitle=sprintf("%s - Overview Fit parameters",myTitle);
    ff=figure('Name',currTitle,'NumberTitle','off');
    if ( ~exist('iRatio','var') )
        iRatio=0;
    end
    if ( exist('OPpath','var') )
        if ( ~exist('beamPart','var') )
            error("For which particle do you want the OP functions?");
        end
    end
    nRows=size(whatToShow,1);
    nCols=size(whatToShow,2);
    for iRow=1:nRows
        for iCol=1:nCols
            iPlot=(iRow-1)*nCols+iCol;
            subplot(nRows,nCols,iPlot);
            tempSplit=split(whatToShow(iRow,iCol),":");
            kick=tempSplit(1);
            obsVar=tempSplit(2);
            obsLoc=tempSplit(3);
            Xs=FitParams(:,iColX,iRow,iCol); Ys=FitParams(:,iColY,iRow,iCol);
            if ( iRatio == 0 )
                tmpMean=mean(Ys);
                plot(Xs,Ys,'.-', ...
                    [min(Xs) max(Xs)],[tmpMean tmpMean],'k-');
                if (exist('OPpath','var'))
                    hold on;
                    Brhos=FitParams(:,iColXrm,iRow,iCol);
                    OPRMs=getOPRMs(magnetName,whatToShow(iRow,iCol),beamPart,Brhos,OPpath);
                    if ( length(OPRMs)==1 )
                        OPRMs=ones(length(Brhos),1)*OPRMs;
                    end
                    plot(Xs,OPRMs,'.-');
                end
                ylabel(labelY);
                legend('MADX',sprintf('MADX ave= %g',tmpMean),'OP','Location','best');
                fprintf("...MADX average: %s\n",tmpMean);
            else
                if (~exist('OPpath','var'))
                    error("cannot plot ratio wrt OP if you don't pass me OP function...");
                end
                Brhos=FitParams(:,iColXrm,iRow,iCol);
                OPRMs=getOPRMs(magnetName,whatToShow(iRow,iCol),beamPart,Brhos,OPpath);
                Ys=Ys./OPRMs;
                tmpMean=mean(Ys);
                plot(Xs,Ys,'.-', ...
                    [min(Xs) max(Xs)],[tmpMean tmpMean],'k-');
                ylabel(sprintf("R_{MAD/OP} of %s",labelY));
                legend('R',sprintf('R ave= %g',tmpMean),'Location','best');
                fprintf("...average R: %s\n",tmpMean);
            end
            grid on;
            xlabel(labelX);
            title(sprintf("%s - %s - %s",obsVar,kick,LabelMe(obsLoc)));
        end
    end
    sgtitle(currTitle);
end

%% xls file handling

function writeXLStable(FitParams,headerFitParams,xlsFileName,outputPath,whatToShow,iDelete)
    origFolder=cd(outputPath);
    fprintf("writing fitting parameters to %s ...\n",xlsFileName);
    iDeleteUsr=1;
    if ( exist('iDelete','var') )
        iDeleteUsr=iDelete;
    end
    if ( iDeleteUsr~=0 )
        delete(xlsFileName);
    end
    nRows=size(whatToShow,1);
    nCols=size(whatToShow,2);
    for iRow=1:nRows
        for iCol=1:nCols
            sheetName=WTS2SN(whatToShow(iRow,iCol));
            writematrix(headerFitParams,xlsFileName,'Sheet',sheetName);
            writematrix(FitParams(:,:,iRow,iCol),xlsFileName,'Sheet',sheetName,'WriteMode','append');
        end
    end
    cd(origFolder);
end

function [FitParams,headerFitParams]=readXLStable(xlsFileName,outputPath,whatToShow)
    origFolder=cd(outputPath);
    fprintf("reading fitting parameters from file %s ...\n",xlsFileName);
    % parse data
    [~,sheet_names]=xlsfinfo(xlsFileName);
    nSpreadSheets=numel(sheet_names);
    nRowsMax=0;
    nColsMax=0;
    for k=1:nSpreadSheets
        [data{k}, myheader]=xlsread(xlsFileName,sheet_names{k});
        if ( find(whatToShow==SM2WTS(sheet_names{k})) >0 )
            if (size(data{k},1)>nRowsMax)
                nRowsMax=size(data{k},1);
            end
            if (size(data{k},2)>nColsMax)
                nColsMax=size(data{k},2);
                headerFitParams=string(myheader);
            end
        end
    end
    % converting into table
    FitParams=zeros(nRowsMax,nColsMax,size(whatToShow,1),size(whatToShow,2));
    for k=1:nSpreadSheets
        if ( find(whatToShow==SM2WTS(sheet_names{k})) >0 )
            [row, column] = find(whatToShow==SM2WTS(sheet_names{k}));
            FitParams(1:size(data{k},1),1:size(data{k},2),row,column)=data{k};
        end
    end
    cd(origFolder);
end

%% Standards

function [whatToShow,fitOrders]=getStandards(magnetName,mode,combo)
% load OP standards for predicting RMs
% input:
% - magnetName: MADX name;
% - mode: "RM" (Response matrix);
% - combo: "REGULAR" (HKICK->x@observation point,VKICK->y@observation point);
%          "COUPLED" (VKICK->x@observation point,HKICK->y@observation point);
% output:
% - whatToShow: matrix (2x2 or 2x3) with actual combinations of kicks to
%               vary and observables;
% - fitOrders: matrix identical to whatToShow expressing the order of the
%              fitting polynom (1: linera, 2: quadratic, etc...);
    if ( ~exist('mode','var') )
        mode="RM";
    end
    if ( ~exist('combo','var') )
        combo="REGULAR";
    end
    if ( strcmpi(combo,"REGULAR") )
        if ( startsWith(upper(magnetName),"T") )
            whatToShow=[ ...
                "HKICK:X:T2_021B_SFH" "HKICK:X:T2_032A_MOB" ; ...
                "VKICK:Y:T2_021B_SFH" "VKICK:Y:T2_032A_MOB" ; ...
            ];
            fitOrders=[ ...
                1 1 ; ...
                1 1 ; ...
            ];
        elseif ( startsWith(upper(magnetName),"U") )
            whatToShow=[ ...
                "HKICK:X:U2_019B_SFH" "HKICK:X:U2_029A_MOB" ; ...
                "VKICK:Y:U2_019B_SFH" "VKICK:Y:U2_029A_MOB" ; ...
            ];
            fitOrders=[ ...
                1 1 ; ...
                1 1 ; ...
            ];
        elseif ( startsWith(upper(magnetName),"V") )
            whatToShow=[ ...
                "HKICK:X:V2_019B_SFH" "HKICK:X:V2_029B_NZL" "HKICK:X:V2_031A_MOB" ; ...
                "VKICK:Y:V2_019B_SFH" "VKICK:Y:V2_029B_NZL" "VKICK:Y:V2_031A_MOB" ; ...
            ];
            fitOrders=[ ...
                1 1 1 ; ...
                1 1 1 ; ...
            ];
        elseif ( startsWith(upper(magnetName),"Z") )
            whatToShow=[ ...
                "HKICK:X:Z2_021B_SFH" "HKICK:X:Z2_032A_MOB" ; ...
                "VKICK:Y:Z2_021B_SFH" "VKICK:Y:Z2_032A_MOB" ; ...
            ];
            fitOrders=[ ...
                1 1 ; ...
                1 1 ; ...
            ];
        else
            error("Unable to detect line of element %s",magnetName);
        end
    elseif ( strcmpi(combo,"COUPLED") )
        if ( startsWith(upper(magnetName),"T") )
            whatToShow=[ ...
                "VKICK:X:T2_021B_SFH" "VKICK:X:T2_032A_MOB" ; ...
                "HKICK:Y:T2_021B_SFH" "HKICK:Y:T2_032A_MOB" ; ...
            ];
            fitOrders=[ ...
                1 1 ; ...
                1 1 ; ...
            ];
        elseif ( startsWith(upper(magnetName),"U") )
            whatToShow=[ ...
                "VKICK:X:U2_019B_SFH" "VKICK:X:U2_029A_MOB" ; ...
                "HKICK:Y:U2_019B_SFH" "HKICK:Y:U2_029A_MOB" ; ...
            ];
            fitOrders=[ ...
                1 1 ; ...
                1 1 ; ...
            ];
        elseif ( startsWith(upper(magnetName),"V") )
            whatToShow=[ ...
                "VKICK:X:V2_019B_SFH" "VKICK:X:V2_029B_NZL" "VKICK:X:V2_031A_MOB" ; ...
                "HKICK:Y:V2_019B_SFH" "HKICK:Y:V2_029B_NZL" "HKICK:Y:V2_031A_MOB" ; ...
            ];
            fitOrders=[ ...
                1 2 2 ; ...
                1 2 2 ; ...
            ];
        elseif ( startsWith(upper(magnetName),"Z") )
            whatToShow=[ ...
                "VKICK:X:Z2_021B_SFH" "VKICK:X:Z2_032A_MOB" ; ...
                "HKICK:Y:Z2_021B_SFH" "HKICK:Y:Z2_032A_MOB" ; ...
            ];
            fitOrders=[ ...
                1 1 ; ...
                1 1 ; ...
            ];
        else
            error("Unable to detect line of element %s",magnetName);
        end
    else
        error("Non-standard request of observations: %s",combo);
    end
end

function OPRMs=getOPRMs(magnetName,whatToShow,beamPart,Brhos,OPpath)
% load appropriate OP RMs
    OPRMs=zeros(length(Brhos),1);
    if ( startsWith(magnetName,"V") )
        addpath(sprintf("%s\%s",OPpath,"\Area dati MD\00Steering\SteeringPazienti\SteeringHEBTlineaV_4.0\RM")); % Line V
        if ( strcmp(beamPart,"PROTON") )
            RmProtoniLineaVfuoco10
        else
            RmCarbonioLineaVfuoco6
        end
        switch whatToShow
            case "HKICK:X:V2_031A_MOB"
                switch magnetName
                    case "V1_044A_CEB"
                        OPRMs=NORM_ISO_Corr1H(Brhos);
                    case "V2_013A_CEB"
                        OPRMs=NORM_ISO_Corr2H(Brhos);
                    otherwise
                        error("Unable to associate an OP RM to %s!",magnetName);
                end
            case "HKICK:X:V2_029B_NZL"
                warning("No OP RM for %s and %s",magnetName,whatToShow);
            case "HKICK:X:V2_019B_SFH"
                switch magnetName
                    case "V1_044A_CEB"
                        OPRMs=NORM_SFH_Corr1H(Brhos);
                    case "V2_013A_CEB"
                        OPRMs=NORM_SFH_Corr2H(Brhos);
                    otherwise
                        error("Unable to associate an OP RM to %s!",magnetName);
                end
            case "VKICK:Y:V2_031A_MOB"
                switch magnetName
                    case "V1_044A_CEB"
                        OPRMs=NORM_ISO_Corr1V(Brhos);
                    case "V2_013A_CEB"
                        OPRMs=NORM_ISO_Corr2V(Brhos);
                    otherwise
                        error("Unable to associate an OP RM to %s!",magnetName);
                end
            case "VKICK:Y:V2_029B_NZL"
                switch magnetName
                    case "V1_044A_CEB"
                        OPRMs=NORM_NZF_Corr1V(Brhos);
                    case "V2_013A_CEB"
                        OPRMs=NORM_NZF_Corr2V(Brhos);
                    otherwise
                        error("Unable to associate an OP RM to %s!",magnetName);
                end
            case "VKICK:Y:V2_019B_SFH"
                switch magnetName
                    case "V1_044A_CEB"
                        OPRMs=NORM_SFH_Corr1V(Brhos);
                    case "V2_013A_CEB"
                        OPRMs=NORM_SFH_Corr2V(Brhos);
                    otherwise
                        error("Unable to associate an OP RM to %s!",magnetName);
                end
            otherwise
        end
    else
        addpath(sprintf("%s\%s",OPpath,"\Area dati MD\00Steering\SteeringPazienti\SteeringHEBT_4.0\RM")); % Lines T,U,Z
        if ( startsWith(magnetName,"T") )
            if ( strcmp(beamPart,"PROTON") )
                RmProtoniLineaT_FG
            else
                RmCarbonioLineaT_FP
            end
            switch whatToShow
                case "HKICK:X:T2_032A_MOB"
                    switch magnetName
                        case "T1_011A_CEB"
                            OPRMs=ISOCorr1H(Brhos);
                        case "T2_015A_CEB"
                            OPRMs=ISOCorr2H(Brhos);
                        otherwise
                            error("Unable to associate an OP RM to %s!",magnetName);
                    end
                case "HKICK:X:T2_021B_SFH"
                    switch magnetName
                        case "T1_011A_CEB"
                            OPRMs=SFHCorr1H(Brhos);
                        case "T2_015A_CEB"
                            OPRMs=SFHCorr2H(Brhos);
                        otherwise
                            error("Unable to associate an OP RM to %s!",magnetName);
                    end
                case "VKICK:Y:T2_032A_MOB"
                    switch magnetName
                        case "T1_011A_CEB"
                            OPRMs=ISOCorr1V(Brhos);
                        case "T2_015A_CEB"
                            OPRMs=ISOCorr2V(Brhos);
                        otherwise
                            error("Unable to associate an OP RM to %s!",magnetName);
                    end
                case "VKICK:Y:T2_021B_SFH"
                    switch magnetName
                        case "T1_011A_CEB"
                            OPRMs=SFHCorr1V(Brhos);
                        case "T2_015A_CEB"
                            OPRMs=SFHCorr2V(Brhos);
                        otherwise
                            error("Unable to associate an OP RM to %s!",magnetName);
                    end
                otherwise
                    error("Unable to associate an OP RM to %s!",magnetName);
            end
        elseif ( startsWith(magnetName,"U") )
            if ( strcmp(beamPart,"PROTON") )
                RmProtoniLineaU_FG
            else
                RmCarbonioLineaU_FP
            end
            switch whatToShow
                case "HKICK:X:U2_029A_MOB"
                    switch magnetName
                        case "U1_023A_CEB"
                            OPRMs=ISOCorr1H(Brhos);
                        case "U2_013A_CEB"
                            OPRMs=ISOCorr2H(Brhos);
                        otherwise
                            error("Unable to associate an OP RM to %s!",magnetName);
                    end
                case "HKICK:X:U2_019B_SFH"
                    switch magnetName
                        case "U1_023A_CEB"
                            OPRMs=SFHCorr1H(Brhos);
                        case "U2_013A_CEB"
                            OPRMs=SFHCorr2H(Brhos);
                        otherwise
                            error("Unable to associate an OP RM to %s!",magnetName);
                    end
                case "VKICK:Y:U2_029A_MOB"
                    switch magnetName
                        case "U1_023A_CEB"
                            OPRMs=ISOCorr1V(Brhos);
                        case "U2_013A_CEB"
                            OPRMs=ISOCorr2V(Brhos);
                        otherwise
                            error("Unable to associate an OP RM to %s!",magnetName);
                    end
                case "VKICK:Y:U2_019B_SFH"
                    switch magnetName
                        case "U1_023A_CEB"
                            OPRMs=SFHCorr1V(Brhos);
                        case "U2_013A_CEB"
                            OPRMs=SFHCorr2V(Brhos);
                        otherwise
                            error("Unable to associate an OP RM to %s!",magnetName);
                    end
                otherwise
                    error("Unable to associate an OP RM to %s!",magnetName);
            end
        elseif ( startsWith(magnetName,"Z") )
            if ( strcmp(beamPart,"PROTON") )
                RmProtoniLineaZ_FG
            else
                RmCarbonioLineaZ_FP
            end
            switch whatToShow
                case "HKICK:X:Z2_032A_MOB"
                    switch magnetName
                        case "Z1_011A_CEB"
                            OPRMs=ISOCorr1H(Brhos);
                        case "Z2_015A_CEB"
                            OPRMs=ISOCorr2H(Brhos);
                        otherwise
                            error("Unable to associate an OP RM to %s!",magnetName);
                    end
                case "HKICK:X:Z2_021B_SFH"
                    switch magnetName
                        case "Z1_011A_CEB"
                            OPRMs=SFHCorr1H(Brhos);
                        case "Z2_015A_CEB"
                            OPRMs=SFHCorr2H(Brhos);
                        otherwise
                            error("Unable to associate an OP RM to %s!",magnetName);
                    end
                case "VKICK:Y:Z2_032A_MOB"
                    switch magnetName
                        case "Z1_011A_CEB"
                            OPRMs=ISOCorr1V(Brhos);
                        case "Z2_015A_CEB"
                            OPRMs=ISOCorr2V(Brhos);
                        otherwise
                            error("Unable to associate an OP RM to %s!",magnetName);
                    end
                case "VKICK:Y:Z2_021B_SFH"
                    switch magnetName
                        case "Z1_011A_CEB"
                            OPRMs=SFHCorr1V(Brhos);
                        case "Z2_015A_CEB"
                            OPRMs=SFHCorr2V(Brhos);
                        otherwise
                            error("Unable to associate an OP RM to %s!",magnetName);
                    end
                otherwise
                    error("Unable to associate an OP RM to %s!",magnetName);
            end
        end
    end
end

function [multFact,unit]=FactorMe(VarName)
    multFact=1.0;
    unit="m";
    switch upper(VarName)
        case {"X","Y"}
            multFact=1000.0;   
            unit="mm";
    end
end

function unit=MADXunits(VarName)
    switch upper(VarName)
        case {"X","Y","BETX","BETY","DX","DY"}
            unit="m";
        case {"PX","PY","ALFX","ALFY","DPX","DPY"}
            unit="1";
        case {"MUX","MUY"}
            unit="rad";
        otherwise
            error("unable to identify MADX var %s",VarName);
    end
end

function [properties,observations]=WTS_2_PROP_OBS(whatToShow)
% convert info in whatToShow (i.e. which kicks and observables) into
%   requests to MADX, e.g.
% properties=[ "HKICK" "VKICK" ];
% observations=[ "X:V2_031A_MOB" "Y:V2_031A_MOB" "X:V2_029B_NZL" "Y:V2_029B_NZL" "X:V2_019B_SFH" "Y:V2_019B_SFH" ];
    allSplit=split(whatToShow,":");
    properties=unique(allSplit(:,:,1))';
    observations=strings(1,size(allSplit,1)*size(allSplit,2));
    kk=0;
    for ii=1:size(allSplit,1)
        for jj=1:size(allSplit,2)
            kk=kk+1;
            observations(kk)=sprintf("%s:%s",allSplit(ii,jj,2),allSplit(ii,jj,3));
        end
    end
end

%% Formatting stuff

function [iCol,label]=getLabels(varName,iCols,labels,units)
    iCol=iCols(upper(labels)==upper(varName));
    label=sprintf("%s [%s]",labels(upper(labels)==upper(varName)),units(upper(labels)==upper(varName)));
end

function sheetName=WTS2SN(whatToShow)
    sheetName=strrep(whatToShow,":","__");
end

function whatToShow=SM2WTS(sheetName)
    whatToShow=strrep(sheetName,"__",":");
end

function [xlsFileName,MADXfileNames,RMfileNames]=fileNames(magnetName,properties)
    xlsFileName=sprintf("%s__RM.xlsx",magnetName);
    MADXfileNames=strings(4,1);
    MADXfileNames(1)="hebt.madx";              % template master 
    MADXfileNames(2)="hebt_matlab.madx";       % actual master
    MADXfileNames(3)="compute_RM.cmdx";        % template commands
    MADXfileNames(4)="compute_RM_matlab.cmdx"; % actual commands
    RMfileNames=strings(length(properties),1);
    for iProperty=1:length(properties)
        RMfileNames(iProperty)=sprintf("%s__%s__rm.csv",lower(magnetName),lower(properties(iProperty)));
    end
end
