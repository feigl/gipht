%     Initialize your variables, and call xmlread to obtain the document node:

%      findLabel = 'Plot Tools';
    %findLabel = 'Coordinate1';
    findLabel = 'DATA_TYPE';
    findCbk = '';

%      xDoc = xmlread(fullfile(matlabroot,'toolbox','matlab','general','info.xml'));
   xDoc = xmlread('phsig.cor.geo.xml');

%     Find all the listitem elements. The getElementsByTagName method returns a deep list that contains information about the child nodes:

%     allListitems = xDoc.getElementsByTagName('listitem');
    allListitems = xDoc.getElementsByTagName('property');

%         Note:   Lists returned by DOM methods use zero-based indexing.
% 
%     For each listitem, compare the text for the label element to the text
%     you want to find. When you locate the correct label, get the callback
%     text:

    for k = 0:allListitems.getLength-1
       thisListitem = allListitems.item(k)
       
       % Get the label element. In this file, each
       % listitem contains only one label.
       thisList = thisListitem.getElementsByTagName('label');
       thisElement = thisList.item(0);

%        % Check whether this is the label you want.
%        % The text is in the first child node.
%        if strcmp(thisElement.getFirstChild.getData, findLabel)
%            thisList = thisListitem.getElementsByTagName('callback');
%            thisElement = thisList.item(0);
%            findCbk = char(thisElement.getFirstChild.getData);
%            break;
%        end
       
    end

%     Display the final results:

    if ~isempty(findCbk)
        msg = sprintf('Item "%s" has a callback of "%s."',...
                      findLabel, findCbk);
    else
        msg = sprintf('Did not find the "%s" item.', findLabel);
    end
    disp(msg);

