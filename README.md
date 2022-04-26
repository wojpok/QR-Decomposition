# Analiza Numeryczna (M) - Pracownia 2
## Dekompozycja QR

Rozkład macierzy jest narzędziem umożliwiającym przekształcenie macierzy do postaci iloczynu kilku macierzy o nowych właściwościach. Zagadnienie to ma zastosowanie w numerycznym rozwiązywaniu układów liniowych. Jednym z rozważanych rozkładów jest rozkład QR. Polega on na przedstawieniu macierzy wejściowej w postaci iloczynu macierzy ortonormalnej Q oraz macierzy górnotrójkątnej R. Stosowalność rozkładu QR w rozwiązywaniu układów równań liniowych wynika bezpośrednio z własności Q oraz R. Rozważmy układ
$$M\vec{x} = \vec{y}$$
Stosując rozkład QR
$$QR\vec{x} = \vec{y}$$
$$Q^{-1}QR\vec{x}=Q^{-1}\vec{y} $$
$$R\vec{x}=Q^{-1}\vec{y} $$
Ortonormalność Q oznacza, że zachodzi $Q^TQ=QQ^T=Id$, czyli $Q^{-1}=Q^T$. 
$$R\vec{x}=Q^T\vec{y}$$
Rozwiązanie takiego układu jest już proste, ponieważ po prawej stronie jest wektor wartości, a po lewej macierz górnotrójkątna przemnożona przez wektor niewiadomych. Można to zrobić metodą *podstawiania wstecz*. 
