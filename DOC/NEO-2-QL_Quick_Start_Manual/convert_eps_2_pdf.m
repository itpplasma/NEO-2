function convert_eps_2_pdf(pr_name)
path = getenv('PATH');
if strfind(computer,'MAC')            
    setenv('PATH',[path,':/usr/local/bin:/usr/texbin']);
    system(['epstopdf ',pr_name]);
else
    system(['/usr/bin/epstopdf ',pr_name]);
end
setenv('PATH',path)