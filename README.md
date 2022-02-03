# POA : Pipeline de Otimização de Antígenos (Antigen Optimization Pipeline)


## 1. INTRODUÇÃO

O Pipeline de Otimização de Antígenos (POA) é um pipeline semiautomático em python para análise, organização, recuperação e triagem de epítopos resultantes da predição de peptídeos antigênicos em ferramentas web. 

A partir dos resultados obtidos em ferramentas de predição para células B, células T-citotóxicas e células T-auxiliares, o POA identifica os melhores alvos e cria um banco de dados com estes antígenos. Além disso, o POA faz a identificação dos peptídeos com alta similaridade de resíduos entre espécies em um determinado grupo de organismos, com o auxílio da ferramenta Epitope Conservancy Analysis (http://tools.iedb.org/conservancy/), e a caracterização da localização dos epítopos em porções transmembranares, com o auxílio da ferramenta pyTMHMM (https://github.com/bosborne/pyTMHMM).


### 1.1 Predição de Peptídeos Antigênicos

O screening de alto rendimento, o estudo de sequências genômicas de organismos patogênicos e a aplicação de técnicas de vacinologia reversa têm possibilitado a construção de vacinas e testes diagnósticos de forma rápida e eficiente. Neste contexto, o uso de abordagens in silico que integram diversas técnicas e conhecimentos, como a construção de multi-epitopos, tem se tornando cada vez mais presente na produção de produtos biotecnológicos e intervenções imuno-terapêuticas contra muitas doenças virais, bacterianas e parasitárias ([Dar et al., 2021](https://www.nature.com/articles/s41598-021-90868-2); [Enayatkhani et al., 2021](https://www.tandfonline.com/doi/full/10.1080/07391102.2020.1756411?casa_token=Ctu83G9J4aEAAAAA%3Aocl_27U5qwcFn-WXSFZIaaSmChfq-hC5ggb281Z1hGuocFuOXH65acPngTeO6HavrexguLxkdEbfAVQ); [Shey et al., 2019](https://www.nature.com/articles/s41598-019-40833-x)).


Com o POA você organiza os resultados das predições in silico de epítopos de proteínas imunogênicas, obtém uma melhor visualização destes resultados (podendo ser um enorme volume de dados) e otimiza o tempo de análise destas predições. Além disso, direciona a análise fazendo a triagem de epítopos únicos (ou compartilhados) num grupo de organismos e epítopos expostos na estrutura membranar.


Para epítopos de células B, o POA possui suporte para análise do resultado de predição da ferramenta [Bepipred 2.0](https://services.healthtech.dtu.dk/service.php?BepiPred-2.0) e [Predicted Antigenic Peptides/IMED](http://imed.med.ucm.es/Tools/antigenic.pl).
Para epítopos de células T-citotóxica, o POA possui suporte para análise do resultado de predição da ferramenta [NetCTL - 1.2](https://services.healthtech.dtu.dk/service.php?NetCTL-1.2).
Para epítopos de células T-auxiliares, o POA possui suporte para análise do resultado de predição da ferramenta [MHC-II Binding Predictions](http://tools.iedb.org/mhcii/).
Resultados de outros preditores podem também ser analisados pelo POA quando previamente organizados em arquivos fasta padronizados (ver seção 3.1.1 - Arquivos de entrada no POA1).



## 2. INSTALAÇÃO

No momento há apenas uma maneira de instalar o POA v1.0, instalando-o manualmente. Para que ele funcione é necessário que todas as dependências também sejam instaladas (instalação manual).
As dependências necessárias estão listadas abaixo:

Python version > 3.6.0

Pandas version 1.3.4

Numpy version 1.21.4

Biopython version 1.79

pyTMHMM


Depois que as dependências forem instaladas, instale o POA v1.0 (Linux) com: git clone https://github.com/UbiratanBatista/POA_Project.git



## 3. USO DO POA


Por se tratar de um pipeline semi-automático, algumas das etapas ainda precisarão ser feitas de forma manual pelo usuário.


A primeira parte do pipeline (POA1) fará a análise dos resultados de predição e criará um banco de dados de antígenos, com base nos valores de rankeamento estabelecidos por cada ferramenta de predição. Ao final, o POA1 organiza os resultados em múltiplos arquivos fasta para serem utilizados na etapa seguinte.


Os arquivos resultantes da primeira fase do POA (POA1) deverão ser submetidos na plataforma web da ferramenta Epitope Conservancy Analysis (http://tools.iedb.org/conservancy/). Os resultados da análise de conservação deverão ser organizados em uma pasta e esta será submetida a segunda parte do pipeline (POA2).


O POA2 analisará os epítopos de maior ou menor conservação (a depender do objetivo da análise) e aplicará sobre estes epítopos a análise de predição de estruturas transmembranares, TMHMM. Por fim retornará os epítopos resultantes da análise de conservação e sua caracterização quanto à posição na membrana (externa, transmembranar e interna).


### 3.1 POA 1


O POA1 recebe os resultados da predição das ferramentas web em diferentes formatos e utiliza um algoritmo que seleciona os peptídeos de interesse e os organiza em matrizes padronizadas. 


#### 3.1.1 Arquivos de entrada (Arquivos dos resultados de predição)


Como pipeline semi-automático, a primeira coisa a se fazer é realizar a predição nas ferramentas web de predição de epítopos (ver seção 1.1). Para isto deve-se submeter via web as sequências das proteínas de interesse nestas ferramentas de predição e seguir os protocolos e parâmetros definidos para cada método e objetivo. O POA1 recebe os diferentes formatos de arquivos de saída de cada método e executa a análise, indicando e organizando os peptídeos preditamente imunogênicos.


###### 3.1.1.1 Bepipred


O arquivo fasta contendo todas as proteínas para a análise deve ser submetido à ferramenta Bepipred. O resultado esperado é algo como isto:

![Captura de tela de 2022-01-26 17-45-11](https://user-images.githubusercontent.com/72517648/151360191-20879962-26e7-4df5-a79e-1f9d79082c8d.png)


O resultado da análise do Bepipred 2.0 é uma página web (html) contendo as proteínas e as regiões preditamente antigênicas (marcadas em E) em cada sequência. Na página de resultados tem a opção de fazer o Download deste resultado em arquivo .json (JSON Summary). O arquivo json contendo toda esta informação será utilizado pelo POA1.


###### 3.1.1.2 Predicting Antigenic Peptides/IMED


Diferentemente do Bepipred, a ferramenta do IMED não recebe mais que uma proteína por vez. Neste caso, as sequências das proteínas deverão ser submetidas uma por vez, definindo sempre o cabeçalho. O resultado da predição é algo como: 

![Captura de tela de 2022-01-26 17-56-13](https://user-images.githubusercontent.com/72517648/151360188-21f130dd-ae97-4067-b7c9-8f138bff95bf.png)


Esta tabela deve ser copiada (ctrl+c) e colada (ctrl+v) em um arquivo de texto (.txt), seguindo o seguinte formato:

![Captura de tela de 2022-01-26 17-58-37](https://user-images.githubusercontent.com/72517648/151360185-4a3b4863-5eab-4024-99d9-5105ac705e53.png)


O arquivo .txt formatado e contendo os epítopos preditos pelo Predicting Antigenic Peptides/IMED será utilizado pelo POA1.


###### 3.1.1.3 NetCTL


O arquivo fasta contendo todas as proteínas para a análise deve ser submetido à ferramenta NetCTL 1.2. O resultado esperado é algo como isto:

![Captura de tela de 2022-01-26 18-02-17](https://user-images.githubusercontent.com/72517648/151360184-6091362a-bcb6-455f-b363-64e575096806.png)


O resultado de saída da análise do NetCTL 1.2 é uma página web (html) contendo os peptídeos resultantes da análise e aqueles preditamente antigênicos (indicados com < - E). Esta página de resultados deverá ser salva em formato .html, O download pode ser feito com a opção “salvar como” através do click do botão direito do mouse sobre a página. ATENÇÃO: Aguardar a página carregar toda a informação para depois efetuar o Download. O arquivo .html contendo toda esta informação será utilizado pelo POA1.


###### 3.1.1.4 MHC-II Binding Predictions


A ferramenta MHC-II Binding Predictions não avalia separadamente as proteínas quando submetidas ao servidor em um único arquivo .fasta. As proteínas deverão ser submetidas separadamente para a análise. O resultado será algo parecido com isto:

![Captura de tela de 2022-01-26 18-07-50](https://user-images.githubusercontent.com/72517648/151360180-72b31616-f40f-4bce-a78c-25124cffc675.png)


Assim como para o NetCTL, o resultado da predição do MHC-II Binding Predictions é uma página .html (escolha a opção “Text file” na seleção “Output format”) contendo os peptídeos, os scores definidos pelos algoritmos usados no servidor e algumas outras informações. Esta página html deverá ser baixada através da opção “salvar como” (click do botão direito do mouse sobre a página). O usuário deverá nomear o arquivo .html resultante da seguinte forma:

<proteína>_<organismo>.html

ATENÇÃO: Esperar a página de resultados carregar por completo para fazer o download.
Como a predição ocorre com uma proteína por vez, todos os arquivos .html deverão ser organizados em uma única pasta (ou diretório). O POA1 utilizará o caminho (path) desta pasta para encontrar os arquivos.


###### 3.1.1.5 Outros preditores


Esta seção é aberta aos resultados de predição de qualquer outro preditor de peptídeos antigênicos, além dos que foram mencionados acima. Os resultados destas outras predições (não se restringe a um único método) devem ser organizados em um arquivo formato fasta. Antes de ser submetido ao POA1, o arquivo fasta contendo os epítopos preditos deverá ser organizado no seguinte formato:

```
><proteína>_<organismo>_<método_preditor>_<ID_NCBI(se_houver)>_<posição_inicial>_<posição_final>

Epitopo1
 
><proteína>_<organismo>_<método_preditor>_<ID_NCBI(se_houver)>_<posição_inicial>_<posição_final>
 
Epitopo2
…
```

#### 3.1.2 Argumentos obrigatórios

 ![3d](https://user-images.githubusercontent.com/72517648/151359262-9553ca42-63ad-4c9e-976f-0ac6704f49df.png)


 __parâmetro -b:__ Arquivo .Json contendo o resultado da predição do Bepipred 2.0 (seção 3.1.1.1);
 
 __parâmetro -p:__ Arquivo .fasta contendo os resultados da predição do PAP/IMED (seção 3.1.1.2);
 
 __parâmetro -m:__ Diretório/pasta contendo todos os arquivos .html resultados da análise do MHC-II Binding Predictions (seção 3.1.1.4);
 
 __parâmetro -n:__ Arquivo .html contendo o resultado da análise do NetCTL 1.2 (seção 3.1.1.3);
 
 __parâmetro -x:__ Arquivo .fasta contendo o resultado da análise de algum outro preditor escolhido pelo usuário, seguindo formatação pré-definida (ver seção 3.1.1.5);

 
 ATENÇÃO: Deverão ser submetidos os resultados de pelo menos um método de predição de epítopos (ver seção 3.1.1).


 Além dos resultados da predição de epitopos (__argumentos n, m, b, p, e x__), para utilizar o POA1 é necessário submeter ao pipeline  (__parâmetro -f__) um arquivo fasta contendo todas as poliproteínas (ou proteínas, caso tenham sido utilizadas em sequências separadas) usadas para a predição de epítopos. Este arquivo deve seguir a seguinte formatação:


 Para poliproteínas completas:



 ```
 >polyprotein_<nome_do_organismo1>_<identificador1_NCBI(se houver)>...

 Sequencia de aminoácidos1 ...

 >polyprotein_<nome_do_organismo2>_<identificador2_NCBI(se houver)>...

 Sequencia de aminoácidos2 ...
 ...

 ```

 Para as análises realizadas, a informação do cabeçalho será identificada a partir do símbolo separador (“_”).


 Além disso, deverá ser designada obrigatoriamente uma pasta para receber os resultados da análise (__parâmetro -d__).



#### 3.1.3 Argumentos opcionais

  ![3d](https://user-images.githubusercontent.com/72517648/151359262-9553ca42-63ad-4c9e-976f-0ac6704f49df.png)


 __parâmetro -mhla:__ Tipo de alelo de HLA que estará sendo investigado na triagem do resultado da análise de ligantes do MHC-II, alelos humanos HLA-DP, HLA-DQ e HLA-DR (default = DR).
 
 __parâmetro -mic:__ Threshold do IC-50 do algoritmo NN_align 2.3 usado como referência do POA1 para seleção dos epitopos. Valores de IC50 < 50 nM admitem alta afinidade de ligação do epitopo com o MHCII, IC50 < 500 nM é considerado uma ligação de afinidade intermediária e IC50 < 5000 nM caracterizam ligações de baixa afinidade (default = 50).
 
 __parâmetro -bmin:__ Valor mínimo para tamanho de epitopos preditos pela ferramenta Bepipred (default = 0, ou seja não há seleção de tamanho mínimo para os epitopos Bepipred).
 
 __parâmetro -bmax:__ Valor máximo para tamanho dos epitopos preditos pela ferramenta Bepipred (default = 0, ou seja não há seleção de tamanho máximo para os epitopos Bepipred).
 
 __parâmetro -pmin:__ Valor mínimo para tamanho de epitopos preditos pela ferramenta PAP/IMED (default = 0, ou seja não há seleção de tamanho mínimo para os epitopos PAP/IMED).
 
 __parâmetro -pmax:__ Valor máximo para tamanho dos epitopos preditos pela ferramenta PAP/IMED  (default = 0, ou seja não há seleção de tamanho máximo para os epitopos PAP/IMED).
 
 __parâmetro -xmin:__ Valor mínimo para tamanho de epitopos preditos por outros preditores que não aqueles mencionados na seção 3.1.1  (default = 0, ou seja não há seleção de tamanho para estes epítopos).
 
 __parâmetro -xmax:__ Valor máximo para tamanho dos epitopos preditos por outros preditores que não aqueles mencionados na seção 3.1.1 (default = 0, ou seja não há seleção de tamanho para estes epítopos).
 
 __parâmetro -e:__ Opção para que os resultados também sejam organizados em planilha Excel (.xslx), além dos arquivos fastas que são gerados obrigatoriamente (default = n).



#### 3.1.4 Rodando o POA1

Com os arquivos dos resultados das predições em mãos, pode-se executar o POA1.

![2](https://user-images.githubusercontent.com/72517648/151359259-0473b328-cb54-423b-ac3c-f49cefc94213.png)

#### 3.1.5 Arquivos de Saída 

Como resultado da análise do POA1, os epítopos e as proteínas serão organizados em arquivos fasta para serem submetidos a ferramenta web Epitope Analysis Conservancy. O fastas de epítopos serão separados pela espécie do organismo, ou seja, ao final da análise todos os epítopos da espécie X serão organizados em um único fasta. O fasta contendo estes epitopos será nomeado:


X_epitope.fasta (<sp_do_organismo>_epitope.fasta)


Na análise de conservação, quando comparados às sequencias de proteínas da espécie em que foram retirados, estes epítopos darão match de 100% com as proteínas em questão. Desta forma, para evitar redundâncias na comparação entre epitopos-proteínas, o POA1 também cria um fasta contendo todas as proteínas presentes na análise, exceto aquelas do organismo em que os próprios epitopos foram retirados. Este arquivo será nomeado:


X_protein.fasta (<sp_do_organismo_não_incluso_na_análise>_epitope.fasta)


O arquivo em questão não possui as proteínas do organismo X. Este arquivo deverá ser submetido a análise de conservação junto (e somente) com o arquivo que o acompanha, o arquivo de epitopos X_epitope.fasta.


O POA1 constrói ainda um relatório de análises contendo algumas estatísticas dos epitopos analisados. Contêm algumas informações como número de proteínas, número de epitopos e quantidade maior e menor de epitopos dentre as proteínas. Além disso, o pipeline permite a exportação dos epitopos em formato .xslx (Excel) (ver seção 3.1.3 - argumentos opcionais).


### 3.2 POA2


Após a análise de conservação feita pela ferramenta Epitope Conservancy Analysis (IEDB) com os arquivos resultantes do POA1, os dados de conservação para cada peptídeo já podem ser trabalhados. O POA2 recebe os resultados da predição de análise de conservação e seleciona os peptídeos com base em valores de threshold estabelecidos. Além disso, utiliza a ferramenta pyTMHMM  para classificar o resultado quanto à topologia de membrana, a partir da predição de hélices transmembranares. Os peptídeos caracterizados são armazenados em arquivos de planilha (.xlsx) e em arquivos fasta.


#### 3.2.1 Arquivos de entrada


Como discutido, os resultados da primeira metade do pipeline (POA1) devem ser submetidos a análise de conservação na ferramenta online Epitope Conservancy Analysis (IEDB) (http://tools.iedb.org/conservancy/). A análise de conservação deverá ser feita para todos os fastas de epítopos (separados por espécie do organismo no POA1). Os resultados são tabelas (formato .csv) que deverão ser baixadas e armazenadas em um único diretório. O POA2 recebe e avalia cada um dos arquivos fasta presentes no diretório, seleciona os epitopos com base no valor de conservação definido pelo usuário e organiza em uma planilha os epítopos que são mais ou menos conservados (ver seção 3.2.2) e sua classificação quanto a topologia da membrana. 


#### 3.2.2 Argumentos Obrigatórios

 ![Design sem nome](https://user-images.githubusercontent.com/72517648/151359269-a50cbc97-443b-4dc4-9082-b1b52527a566.png)


Para utilizar o POA2 será necessário organizar todos os resultados da análise de conservação em um único diretório (__parâmetro -d__). Além disso, o usuário deverá informar qual o threshold de identidade de sequência (Sequence identity threshold) foi utilizado para fazer a predição (__parâmetro -t__) e se esta foi feita visando encontrar epítopos conservados (>=) ou epítopos únicos (<), __parâmetro -g ou -l__, respectivamente. Deverá ser submetido também um arquivo fasta que contenha todas as polyproteínas usadas para a predição de epítopos (__parâmetro -f__). Este arquivo poderá ser igual àquele submetido no POA1, seguindo a mesma formatação (ver seção 3.1.2).


#### 3.2.3 Argumentos Opcionais

![Design sem nome](https://user-images.githubusercontent.com/72517648/151359269-a50cbc97-443b-4dc4-9082-b1b52527a566.png)


 __parâmetro -r:__ Diretório onde os resultados da análise do POA2 (planilhas e ou arquivos fasta) serão organizados;
 
 __parâmetro -imin:__ Threshold mínimo de identidade aceito no resultado da análise de conservação entre os peptídeos e as proteínas avaliadas (default = 60).
 
 __parâmetro -imax:__ Threshold máximo de identidade aceito no resultado da análise de conservação entre os peptídeos e as proteínas avaliadas (default = 100).
 
 __parâmetro -m:__ porcentagem das sequências que tiveram match de identidade acima do threshold de identidade de sequência avaliado (parâmtero -t) (default = 60).
 
 __parâmetro -rf:__ Opções de arquivo fasta contendo os resultados da análise de conservação (default = None). 
 [0] Todos os epitopos resultantes da análise de conservação, sem considerar a classificação por topologia de membrana.
 [1] Os epitopos resultantes da análise de conservação que estiverem em porções expostas da membrana.
 [2] Os epitopos resultantes da análise de conservação que estiverem em porções transmembranares.
 [3] Os epitopos resultantes da análise de conservação que estiverem em porções internas à membrana.


#### 3.2.4 Rodando o POA2


Com os arquivos dos resultados das predições em mãos, pode-se executar o pipeline.
Ex: Análise de conservação para epítopos conservados (-g True)

![dfd](https://user-images.githubusercontent.com/72517648/151359274-2b88b315-0f56-4494-84a8-7025be6e21df.png)


#### 3.2.5 Arquivos de Saída


Como resultado da análise do POA2, todos os epítopos selecionados na análise de conservação são organizados em um único arquivo excel (.xlsx). Além de conter as informações já disponibilizadas na planilha de resultado da ferramenta do IEDB, a tabela conta com colunas que identificam a porções do epitopo quanto a topologia da membrana, em porções externas, transmembranares e internas (representação da porcentagem, valores de 0 a 1). Esta análise de topologia da membrana é feita usando um algoritmo derivado do TMHMM (Sonnhammer et al., 1998), o pyTMHMM. 
Os epitopos triados podem também ser organizados em arquivos fasta, como demonstrado na seção 3.2.3, em __parâmetro -rf__.
