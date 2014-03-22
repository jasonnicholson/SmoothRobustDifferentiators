publish('index.m','outputDir','..\');
publish('robustDiffDocumentation.m','outputDir','..\');
publish('robustDiffOneSideDocumentation.m','outputDir','..\');

directory = which('generateDocumentation.m');
documentationDirectory = regexprep(directory,['\' filesep 'documentation generators' '\' filesep 'generateDocumentation.m'],'');
builddocsearchdb(documentationDirectory);