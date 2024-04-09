//
// Created by Kamilla Faizullina on 08.04.2024.
//

#include "MonteCarlo.h"


MC_Interacting_SAW_XY::MC_Interacting_SAW_XY(long length,
                                             double Probability_Local_Update,
                                             double Probability_Reconnect) {
    p_for_local_update = Probability_Local_Update;
    p_for_reconnect = p_for_local_update + Probability_Reconnect;
    model = new XY_SAW_LongInteraction(length);
}


void MC_Interacting_SAW_XY::run_simulation() {

}