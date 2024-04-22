#include "global.h"

void EnlightenmentScreen_Init(EnlightenmentScreen* this) {
    this->dListStart = NULL;
    this->dListPtr = NULL;
}

void EnlightenmentScreen_Open(EnlightenmentScreen* this, Gfx* dList) {
    this->dListStart = dList;
    this->dListPtr = this->dListStart;
}

void EnlightenmentScreen_Render(EnlightenmentScreen* this) {
    gDPPipeSync(this->dListPtr++);
    gDPSetRenderMode(this->dListPtr++, G_RM_XLU_SURF, G_RM_XLU_SURF2);
    gDPSetCombineMode(this->dListPtr++, G_CC_SHADE, G_CC_SHADE);
    gDPSetEnvColor(this->dListPtr++, 128, 128, 128, 128);  // Gray color, 50% transparency

    // Assuming a standard N64 screen resolution for the overlay box
    gDPFillRectangle(this->dListPtr++, 0, 0, 320, 240);

    gDPPipeSync(this->dListPtr++);
}

Gfx* EnlightenmentScreen_Close(EnlightenmentScreen* this) {
    Gfx* ret = this->dListPtr;
    this->dListStart = NULL;
    this->dListPtr = NULL;
    return ret;
}
