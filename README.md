# Ready_To_Edit_Figures
목적: Adobe에서 간단하게 최종 figure 만들 수 있는 단순화된 코드 제작하기 (Pathway Enrichment, Sample QC, BED peak)
- Input file 넣으면 output file, figure, table까지 다 만들 수 있게.
  - 단, table은 xlsx로 저장하고, figure 만들었을 때 axis text가 하나의 text로 연결되도록 옵션을 설정할 것.
- Pathway enrichment -> GOBP, Reactome, KEGG 등의 결과가 하나의 plot 안에 다 들어가게 하면 됨.
- 추가로, sample overview (QC EDA 등)가 하나의 supple figure로 정리되는 multipanel figure를 제작하기.
- BED file에서 peak 그림을 만들어내는 plot도 레시피를 만들어두기 (Guitarplot 등?)

### Random thoughts.
- Pathway term이 되게 방대해질텐데, 관심있는 term만으로 plot을 만들 수 있도록 subset하려면 어떻게 해야할까? 가령 grepl이나..?
