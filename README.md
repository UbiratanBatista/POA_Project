# POA_Project
Optimization pipeline of predicted putative antigens in web prediction tools. 

# **POA : Pipeline de Otimização de Antígenos (Antigen Optimization Pipeline) 

# Introdução
O Pipeline de Otimização de Antígenos (POA) é um pipeline semiautomático em python para analise, organização, recuperação e triagem de epítopos resultados da predição de peptídeos antigênicos em ferramentas web. 
A partir dos resultados descritos em algumas ferramentas de predição para células B, células T-citotóxicas e células T-auxiliares, o POA identifica os melhores alvos e cria um banco de dados com estes antígenos. O POA faz ainda:
Identificação de peptídeos com alta similaridade de resíduos entre espécies em um grupo de organismos, com o auxílio da ferramenta Epitope Conservancy Analysis.
Caracterização da localização do epítopo em porções transmembranares, com o auxilio da ferramenta TMHMM.

Predição de Peptídeos Antigênicos 
O screening de alto rendimento, o estudo de sequências genômicas de organismos patogênicos e a aplicação de técnicas de vacinologia reversa tem possibilitado a construção de vacinas e testes diagnósticos de forma rápida e eficiente. Neste contexto, abordagens integrativas, como os multi-epitopos, para este fim tem se tornando cada vez mais presente na intervenção farmacológica (??) para muitas doenças virais, bacterianas e parasitárias. Estudos recentes têm corroborado a eficiência de técnicas in sílico para a construção destes insumos biotecnológicos (BUSCAR REFERENCIA).
Com o POA você organiza os resultados das predições in sílico de epítopos de proteínas imunogênicas, obtém uma melhor visualização destes resultados (podendo ser um enorme volume de dados) e otimiza o tempo de análise destas predições. Além disso, direciona a análise fazendo a triagem de epítopos únicos (ou compartilhados) num grupo de organismos e epítopos expostos na estrutura membranar.
Para epítopos de células B, o POA possui suporte para análise do resultado de predição da ferramenta Bepipred 2.0 (https://services.healthtech.dtu.dk/service.php?BepiPred-2.0) e Predicted Antigenic Peptides/IMED (http://imed.med.ucm.es/Tools/antigenic.pl).
Para epítopos de células T-citotóxica, o POA possui suporte para análise do resultado de predição da ferramenta NetCTL - 1.2 (https://services.healthtech.dtu.dk/service.php?NetCTL-1.2).
Para epítopos de células T-auxiliares, o POA possui suporte para análise do resultado de predição da ferramenta MHC-II Binding Predictions (http://tools.iedb.org/mhcii/).
Resultados de outros preditores podem também ser analisados pelo POA quando previamente organizados em fastas padronizados (ver seção arquivos de entrada no POA1).

Instalação
No momento há apenas uma maneira de instalar o POA v1.0, instalando-o manualmente. Para que ele funcione é necessário que todas as dependências também sejam instaladas (instalação manual).
As dependências necessárias estão listadas abaixo:
Python version > 3.6.0
Pandas version 1.3.4
Numpy version 1.21.4
Biopython version 1.79
pyTMHMM

Depois que as dependências forem instaladas, instale o POA v1.0 (Linux) com: 
git clone https://github.com/UbiratanBatista/POA_Project.git

Uso do POA
Por se tratar de um pipeline semi-automático, algumas das etapas ainda precisarão ser feitas de forma manual pelo usuário.
A primeira parte do pipeline (POA1) fará a análise dos resultados de predição e criará um banco de dados de antígenos, com base nos valores estabelecidos por cada ferramenta de predição. Ao final, o POA1 organiza os resultados em múltiplos arquivos fasta para serem utilizados na etapa seguinte.
Os arquivos resultantes da primeira fase do POA (POA1) deverão ser submetidos na plataforma web da ferramenta Epitope Conservancy Analysis (http://tools.iedb.org/conservancy/). Os resultados da análise de conservação deverão ser organizados em uma pasta e esta será submetida a segunda parte do pipeline (POA2).
O POA2 analisará os epítopos de maior ou menor conservação (a depender do objetivo da análise) e aplicará sobre estes epítopos a análise de predição de estruturas transmembranares, TMHMM. Por fim retornará os epítopos resultantes da análise conservação e sua caracterização quanto à posição na membrana (externa, transmembranar e interna).

POA1

O POA1 recebe os resultados da predição das ferramentas web em diferentes formatos e utiliza um algoritmo que seleciona os peptídeos de interesse em os organiza em tabelas padronizadas. 
3.1.1 Argumentos obrigatórios
<imagem dos argumentos obrigatórios>
-Para utilizar o POA1 é necessário organizar todas as proteínas usadas para a predição de epítopos em um único arquivo fasta, parâmetro -f, seguindo a formatação:
><nome_da_proteína1>_<nome_do_organismo1>_<identificador1_NCBI (se houver)>...
Sequencia de aminoácidos ...
><nome_da_proteína2>_<nome_do_organismo2>_<identificador2_NCBI (se houver)>...
Sequencia de aminoácidos ...
...
Para as análises realizadas, a informação do cabeçalho será identificada a partir do símbolo separador (“_”).
- Deverá ser designada obrigatoriamente uma pasta para receber os resultados da análise, parâmetro -d.
-Deverá ser submetido pelo menos um arquivo resultado da predição de epítopos (ver seção 1.1).
3.1.2 Arquivos de entrada (Arquivos dos resultados de predição)
Como pipeline semi-automático, a primeira coisa a se fazer é realizar a predição nas ferramentas web de predição de epítopos (ver seção 1.1). Para isto deve-se submeter as sequencias das proteínas de interesse nas ferramentas de predição (pelo menos uma). O POA recebe os diferentes formatos de arquivos de saída e executa a análise, indicando quais são os peptídeos preditamente imunogênicos.
3.1.2.1 Bepipred
O arquivo fasta contendo todas as proteínas para a análise deve ser submetido à ferramenta Bepipred. O resultado esperado é algo como isto:
<imagem do output Bepipred>
O resultado de saída da análise do Bepipred 2.0 é uma página web (html) contendo as proteínas e as regiões preditamente antigênicas (marcadas em B). Na página de resultados tem a opção de fazer o Download deste resultado em arquivo .json. O arquivo json contendo toda esta informação será utilizado pelo POA.

3.1.2.2 Predicting Antigenic Peptides/IMED
Diferentemente do Bepipred, a ferramenta do IMED não recebe mais que uma proteína por vez. Neste caso, as sequencias das proteínas deverão ser submetidas uma por vez, definindo sempre o cabeçalho. O resultado da predição é algo como: 
<imagem do output PAP/IMED >
Esta tabela deve ser copiada (ctrl+c) e colada (ctrl+v) em um arquivo de texto (.txt), seguindo o seguinte formato:
<imagem do arquivo txt>
O arquivo .txt formatado e contendo os epítopos preditos pelo Predicting Antigenic Peptides/IMED será utilizado pelo POA.

3.1.2.3 NetCTL
O arquivo fasta contendo todas as proteínas para a análise deve ser submetido à ferramenta NetCTL 1.2. O resultado esperado é algo como isto:
<imagem do output NetCTL >
O resultado de saída da análise do NetCTL 1.2 é uma página web (html) contendo os peptídeos resultantes da análise e aqueles preditamente antigênicos (indicados com < - E). Esta página de resultados deverá ser salva em formato .html, O download pode ser feito com a opção salvar como através do click do botão direito do mouse sobre a página. O arquivo html contendo toda esta informação será utilizado pelo POA.

3.1.2.4 MHC-II Binding Predictions
A ferramenta não avalia separadamente as proteínas submetidas para predição se unidas em um único fasta. As proteínas deverão ser submetidas separadamente para a análise do MHC-II Binding Predictions. O resultado será algo parecido com isto:
<imagem do output MHCII>
Assim como para o NetCTL, o resultado da predição do MHC-II Binding Predictions é uma página .html contendo os peptídeos e os scores definidos pelo servidor. Esta página html deverá ser baixada através da opção “salvar como” (click do botão direito do mouse sobre a página). O usuário deverá deverá nomear o arquivo .html resultante da seguinte forma:
<proteína>_<organismo>.html
Como a predição ocorre uma proteína por vez, todos os arquivos html deverão ser organizados em uma única pasta (ou diretório).

3.1.2.5 Outros preditores
Esta seção é aberta aos resultados de predição de qualquer outro preditor de peptídeos antigênicos, além dos que foram mencionados acima. Os resultados destas outras predições (não se restringe a um único método) devem ser organizados em um arquivo formato fasta. O fasta contendo os epítopos preditos deverá ser organizado no seguinte formato:
><proteína>_<organismo>_<método_preditor>_<ID_NCBI(se_houver)>_<posição_inicial>_<posição_final>
Epitopo1
><proteína>_<organismo>_<método_preditor>_<ID_NCBI(se_houver)>_<posição_inicial>_<posição_final>
Epitopo2
…

Argumentos opcionais
Explicar que tem como modificar o output modificando alguns parâmetros
Pontuar em cada método


3.1.3 Rodando o POA1
Com os arquivos dos resultados das predições em mãos, pode-se executar o pipeline.
<print da linha de comando>

3.1.4 Arquivos de Saída 
Como resultado da análise do POA1, os epítopos e as proteínas serão organizados em arquivos fasta para serem submetidos a ferramenta web Epitope Analysis Conservancy. O fastas de epítopos serão separados pela espécie do organismo, ou seja, ao final da análise todos os epítopos da espécie X serão organizados em um único fasta. O fasta contendo estes epitopos será nomeado:
X_epitope.fasta (<sp_do_organismo>_epitope.fasta)
Na análise de conservação, quando comparados às sequencias de proteínas da espécie em que foram retirados, estes epitopos darão match de 100% com as proteínas em questão. Desta forma, para evitar redundâncias na comparação entre epitopos-proteinas, o POA1 também cria um fasta contendo todas as proteínas presentes na ánalise, exceto aquelas do organismo em que os próprios epitopos foram retirados. Este arquivo será nomeado:
X_protein.fasta (<sp_do_organismo_não_incluso_na_análise>_epitope.fasta)
O arquivo em questão não possui as proteínas do organismo X. Este arquivo deverá ser submetido a análise de conservação junto (e somente) com o arquivo que o acompanha, o arquivo de epitopos X_epitope.fasta.
O POA1 constrói ainda um relatório de análises contendo algumas estatísticas dos epitopos analisados. Contêm algumas informações como número de proteínas, número de epitopos e quantidade maior e menor de epitopos dentre as proteínas. Além disso, o pipeline permite a exportação dos epitopos em formato .xslx (Excel) (ver seção argumentos opcionais).


>>> POA2
