# missing-data-pks

Der Code ist in mehreren "Ebenen" organisiert:
1. R Funktionen in dem Ordner "funktionen" und C++ Funktionen im Ordner "funktionen_cpp"
   - "funktionen":
      - genKS.R: simuliert k�nstliche Wissensstruktur
      - simMissingsBlim.R: simuliert basierend auf einer WS samt beta, eta, mu_q und mu_q_ Antwortmuster

   - "funktionen_cpp":
      - BLIM_combi.R: Wrapper Funktion, ruft die C++ Funktionen EMblim.cpp, EMimblim.cpp und EMmissblim.cpp auf, welche entsprechend die EM-Modellsch�tzung durchf�hren.

2. "�bergeordnete" Funktionen
   - data_sim.R: Simuliert Datens�tze basierend auf der simulierten WS mit entsprechenden mu_q und mu_q_ Parametern, ruft simMissingsBlim.R auf.
   - data_sim_emp.R: �quivalent nur hier f�r die empirischen WS mit anderen mu_q und mu_q_ Parametern.
   - fit_models.R: Passt alle Modelle an die Daten an; ruft BLIM_combi.R auf.

3. Simulationsablauf ("flow"): WS simulieren (genKS.R)/einlesen, Daten simulieren (data_sim.R bzw. data_sim_emp.R) und Modelle anpassen (fit_models.R)
   - simWS_flow.R: Ablauf f�r simulierte, k�nstliche WS
   - size_flow.R: Ablauf f�r verschiedene Stichproben mit k�nstlicher WS
   - empWS_flow.R: Ablauf f�r empirische WS (vgl. Ordner fpi_WS)
