#ifndef _Z_BG_JYA_BLOCK_H_
#define _Z_BG_JYA_BLOCK_H_

#include <ultra64.h>
#include <global.h>

typedef struct {
    /* 0x0000 */ Actor actor;
    /* 0x014C */ char unk_14C[0x18];
} BgJyaBlock; // size = 0x0164

extern const ActorInit Bg_Jya_Block_InitVars;

#endif
