# missing-data-pks

Der Code ist in mehreren "Ebenen" organisiert:
1. R Funktionen in dem Ordner "funktionen" und C++ Funktionen im Ordner "funktionen_cpp"
   - "funktionen":
      - genKS.R: simuliert künstliche Wissensstruktur
      - simMissingsBlim.R: simuliert basierend auf einer WS samt beta, eta, mu_q und mu_q_ Antwortmuster

   - "funktionen_cpp":
      - BLIM_combi.R: Wrapper Funktion, ruft die C++ Funktionen EMblim.cpp, EMimblim.cpp und EMmissblim.cpp auf, welche entsprechend die EM-Modellschätzung durchführen.

2. "Übergeordnete" Funktionen
   - data_sim.R: Simuliert Datensätze basierend auf der simulierten WS mit entsprechenden mu_q und mu_q_ Parametern, ruft simMissingsBlim.R auf.
   - data_sim_emp.R: äquivalent nur hier für die empirischen WS mit anderen mu_q und mu_q_ Parametern.
   - fit_models.R: Passt alle Modelle an die Daten an; ruft BLIM_combi.R auf.

3. Simulationsablauf ("flow"): WS simulieren (genKS.R)/einlesen, Daten simulieren (data_sim.R bzw. data_sim_emp.R) und Modelle anpassen (fit_models.R)
   - simWS_flow.R: Ablauf für simulierte, künstliche WS
   - size_flow.R: Ablauf für verschiedene Stichproben mit künstlicher WS
   - empWS_flow.R: Ablauf für empirische WS (vgl. Ordner fpi_WS)
