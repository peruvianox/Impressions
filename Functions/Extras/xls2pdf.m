function [status,message] = xls2pdf(xlFile,sheetNameOrNumber,pdfFile,XLRange)
    %% XLS2PDF. Save Microsoft Excel file into PDF programatically.
    %
    %  [status,message] = xls2pdf(xlFile,sheetNameOrNumber,pdfFile);
    %  [status,message] = xls2pdf(xlFile,sheetNameOrNumber,pdfFile,XLRange);
    %  Any file that can be opened with MS Excel can be saved into PDF.
    %
    %  FILE NAMES MUST BE FULL PATH NAMES i.e.: C:\Users\Desktop\Daniel\example.xlxs
    %  
    %  status is an output variable indicating successful conversion (true) or not (false).
    %  message appears if status is false, indicating the cause of the error.
    %
    %  xlFile is the Microsoft Excel file to be converted, specifying full path
    %  sheetNameOrNumber is the sheet Name (sting) or Number (integer) that is desired to be converted into pdf
    %  pdfFile is the PDF file into which spreadsheet should be converted, specifying full path
    %  XLRange is the range, specified in excel format (e.g. 'A1:B50') to be exported to the PDF file
    %
    %  Example:
    %
    %   [OK,msg] = xls2pdf([pwd filesep 'myExcel.xlsx'],'Sheet1',fullfile(pwd,'PDFs/myExcelConverted.pdf'));
    %   [OK,msg] = xls2pdf([pwd filesep 'myExcel.xlsx'],1,fullfile(pwd,'PDFs/myExcelConverted.pdf'),'A1:B50');
    %
    %
    % To convert form Word to PDF, use the doc2pdf tool:
    %   http://www.mathworks.com/matlabcentral/fileexchange/47308-doc2pdf-m
    %
    % To generate full paths to files, you can use GetFullPath from Matlab
    % Central File Exchange: http://www.mathworks.com/matlabcentral/fileexchange/28249-getfullpath
    
    if nargin~=3 && nargin~=4                                              % If number of arguments is ~=3 && ~=4
       error('Wrong numer of inputs. xls2pdf(xlFile,sheetNameOrNumber,pdfFile,...)'); % Error
    end
    
    status  = false;                              % Initialize status
    message = '';                                 % Initialize message
    
    try
        Excel    = actxserver('Excel.Application');          % Open MS Excel
        Workbook = Excel.Workbooks.Open(xlFile);             % Open file
        Sheet    = Workbook.Sheets.Item(sheetNameOrNumber);  % Select Sheet
        if     nargin==3                                     % If Range is not specified (3 inputs)
            Sheet.ExportAsFixedFormat('xlTypePDF',pdfFile);  % Export Sheet
        else                                                 % If Range is specified (4 inputs)
            Sheet.Range(XLRange).ExportAsFixedFormat('xlTypePDF',pdfFile); % Export Range
        end

        [~, ~, ext] = fileparts(xlFile);                     % Analyze extension
        if strcmpi(ext,'.xlsx')
            Workbook.Close(false)                            % Close workbook withou saving
            Excel.Quit;                                      % Quit Excel
            delete(Excel);                                   % Delete matlab object
        else
            Excel.Quit
            Excel.delete
        end
        status  = true;                                      % Successfully saved
    catch exc
        message = sprintf('%s\n%s',exc.identifier,exc.message); % show message
    end