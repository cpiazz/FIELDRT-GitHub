function FIELDRTwheel(hObject, callbackdata, handles)

tmpSlice = [];

if isempty(callbackdata) % when the callback is empty
    % Just to refresh the axial handler that keeps track of the number of
    % the slice
    handles.currentSlice_ax = handles.currentSlice_ax;
else
    % If the scroll is moved up or down
    if callbackdata.VerticalScrollCount > 0
        tmpSlice = handles.currentSlice_ax + 1;
    elseif callbackdata.VerticalScrollCount < 0
        tmpSlice = handles.currentSlice_ax - 1;
    end
    
    % Some information retrieved
    indexS = handles.planC{end};
    structfieldnum = indexS.structures;
    scanfieldnum = indexS.scan;
    volsize = size(handles.planC{1, scanfieldnum}.scanArray);
    
    % Condition to limit the changing to the image size
    if isempty(tmpSlice)
        
    else
        if tmpSlice < 1 || tmpSlice > size(handles.planC{1, handles.scanfieldnum}.scanArray, 3)
            % Slice number is not changing
        else
            % Slice number is changing
            handles.currentSlice_ax = tmpSlice;
        end
    end
    
    % updating axial view
    axes(handles.axes3)
    
    % Deleting all annotations in the image
    axesHandlesToChildObjects = findobj(gca, 'Type', 'Line');
    if ~isempty(axesHandlesToChildObjects)
        delete(axesHandlesToChildObjects);
    end
    
    % delete(findall(gcf,'type','annotation'))
    
    
    set(handles.axes3.Children, 'cdata', handles.planC{1, handles.scanfieldnum}.scanArray(:, :, handles.currentSlice_ax));
    % set(findobj('Tag','Image'), 'cdata', handles.planC{1, handles.scanfieldnum}.scanArray(:, :, handles.currentSlice_ax));
    % Updating the image and the number of the slice
    set(handles.numslice, 'String', [int2str(handles.currentSlice_ax) '/' int2str(size(handles.planC{1, handles.scanfieldnum}.scanArray, 3))]);
%     set(handles.hShow, 'cdata', handles.planC{1, handles.scanfieldnum}.scanArray(:, :, handles.currentSlice_ax));
    
%     % updating sagittal view
%     axes(handles.axes5)
%     set(handles.axes5.Children, 'cdata', handles.ImgSg(:, :, handles.currentSlice_sag));
%     % Updating the image and the number of the slice
%     set(handles.numslice2, 'String', [int2str(handles.currentSlice_sag) '/' int2str(size(handles.ImgSg, 3))]);
% 
%     % updating coronal view
%     axes(handles.axes6)
%     set(handles.axes6.Children, 'cdata', handles.ImgCr(:, :, handles.currentSlice_cor));
%     % Updating the image and the number of the slice
%     set(handles.numslice3, 'String', [int2str(handles.currentSlice_cor) '/' int2str(size(handles.ImgCr, 3))]);
    
    % imshow(handles.planC{1, scanfieldnum}.scanArray(:, :, handles.currentSlice_ax), []);
    % text = sprintf('Slice number: %d', handles.currentSlice_ax);
    % set(handles.ct_slice_number, 'String', text);
    index_selected = get(handles.listbox1,'Value');
    FIELDRTplotContours(handles, index_selected)
    
    % Updating the text panel with the quantitavie feedback (only Jaccard 2D
    % will be updated)
    % The update is not going to happen when oars or minmax toogle buttons are clicked
    
    if get(handles.oars, 'Value') == 1 || get(handles.accreg, 'Value') == 1
        quantitatfeed = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n',[]);
        
        set(handles.stats, 'String', quantitatfeed);
    else
        % Volumes
        if isempty(index_selected)
            
        else
            volGSvalue = num2str(handles.volumestat_GS{index_selected, 1}.vol);
            
            if size(volGSvalue, 2) < 6 % if the number is too long!
                
            else
                volGSvalue = volGSvalue(1 : 6);
            end
            
            volUSvalue = num2str(handles.volumestat_US{index_selected, 1}.vol);
            
            if size(volUSvalue, 2) < 6 % if the number is too long!
                
            else
                volUSvalue = volUSvalue(1 : 6);
            end
            
            %     % We don't wait this information
            %     % Volume Ratio
            %     volratiovalue = num2str(handles.volume_Ratio{index_selected, 1});
            %
            %     if size(volratiovalue, 2) < 6 % if the number is too long!
            %
            %     else
            %         volratiovalue = volratiovalue(1 : 6);
            %     end
            %
            %     % COM GS
            %     COMGSvalue = num2str([num2str(handles.volumestat_GS{index_selected, 1}.COM(1, 1)) ' ' ...
            %         num2str(handles.volumestat_GS{index_selected, 1}.COM(1, 2)) ' ' ...
            %         num2str(handles.volumestat_GS{index_selected, 1}.COM(1, 3))]);
            %
            %     % COM US
            %     COMUSvalue = num2str([num2str(handles.volumestat_US{index_selected, 1}.COM(1, 1)) ' ' ...
            %         num2str(handles.volumestat_US{index_selected, 1}.COM(1, 2)) ' ' ...
            %         num2str(handles.volumestat_US{index_selected, 1}.COM(1, 3))]);
            
            % 3D Jaccard similarity coefficient
            jaccard3Dvalue = num2str(handles.jaccard3Doutput{index_selected, 1});
            
            if size(jaccard3Dvalue, 2) < 6 % if the number is too long!
                
            else
                jaccard3Dvalue = jaccard3Dvalue(1 : 6);
            end
            
            % 2D Jaccard similarity coefficient
            jaccard2Dvalue = num2str(handles.jaccard2Doutput{index_selected, handles.currentSlice_ax});
            
            if size(jaccard2Dvalue, 2) < 6 % if the number is too long!
                
            else
                jaccard2Dvalue = jaccard2Dvalue(1 : 6);
            end
            
            % Storing all this information in one variable
            quantitatfeed = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n',['Volume GS = ' volGSvalue], ['Volume US = ' volUSvalue], ...
                ['3D Jaccard similarity coefficient = ' jaccard3Dvalue], ...
                ['2D Jaccard similarity coefficient = ' jaccard2Dvalue]);
            %        ['Volume Ratio = ' volratiovalue], ['COM GS = ' COMGSvalue], ['COM US = ' COMUSvalue], ...
            
            set(handles.stats, 'String', quantitatfeed);
        end
    end
    
    
    
    % Updating the number of the current slice
    % display(handles.currentSlice)
    handles.currentSlice_ax = round(handles.currentSlice_ax);
    refreshdata(handles.axes3);
    set(gcf, 'WindowScrollWheelFcn', {@FIELDRTwheel,handles});
    
    guidata(hObject,handles);
end