#include "global.h"

MtxF D_801600C0[60];
char D_80160FC0[0x10];

void func_800A57C0(MtxF* mtx, Struct_800A57C0* arg1, Struct_800A598C* arg2, Vtx* arg3, Vec3f* arg4) {
    Vtx* temp_s2;
    Struct_800A57C0* phi_s1;
    f32 temp_x;
    f32 temp_y;
    f32 temp_z;
    Vec3f sp70;
    Vec3f sp64;

    for (phi_s1 = arg1; phi_s1 < arg1 + arg2->unk_0; phi_s1++) {
        temp_s2 = &arg3[phi_s1->unk_0];
        temp_s2->n.ob[0] = arg4->x;
        temp_s2->n.ob[1] = arg4->y;
        temp_s2->n.ob[2] = arg4->z;
        temp_x = mtx->wx;
        temp_y = mtx->wy;
        temp_z = mtx->wz;
        mtx->wz = 0.0f;
        mtx->wy = 0.0f;
        mtx->wx = 0.0f;
        sp64.x = phi_s1->unk_6;
        sp64.y = phi_s1->unk_7;
        sp64.z = phi_s1->unk_8;
        SkinMatrix_Vec3fMtxFMultXYZ(mtx, &sp64, &sp70);
        temp_s2->n.n[0] = sp70.x;
        temp_s2->n.n[1] = sp70.y;
        temp_s2->n.n[2] = sp70.z;
        mtx->wx = temp_x;
        mtx->wy = temp_y;
        mtx->wz = temp_z;
    }
}

#ifdef NON_MATCHING
// Regalloc, Stack
void func_800A598C(GraphicsContext* gfxCtx, PSkinAwb* skin, s32 limbIndex, s32 arg3) {
    SkinLimb** limbs;
    SkinLimb* limb;
    Struct_800A5E28* data;
    Struct_800A598C* phi_s6;
    s32 temp_a0_2;
    SkinAvb* avb;
    Struct_800A598C* new_var;
    Struct_800A57C0* temp_v0;
    Struct_800A598C_2* temp_s3;
    Vtx* spEC;
    Vec3f spDC;
    Vec3f spD0;

    OPEN_DISPS(gfxCtx, "../z_skin.c", 254);

    limbs = (SkinLimb**)SEGMENTED_TO_VIRTUAL(skin->skeletonHeader->segment);
    limb = (SkinLimb*)SEGMENTED_TO_VIRTUAL(limbs[limbIndex]);
    data = (Struct_800A5E28*)SEGMENTED_TO_VIRTUAL(limb->segment);
    phi_s6 = &((Struct_800A598C*)SEGMENTED_TO_VIRTUAL(data->unk_4))[0];
    avb = &skin->avbTbl[limbIndex];
    spEC = avb->buf[avb->unk_0];
    new_var = &phi_s6[data->unk_2];
    while (phi_s6 < new_var) {
        temp_a0_2 = phi_s6->unk_2;
        temp_v0 = (Struct_800A57C0*)SEGMENTED_TO_VIRTUAL(phi_s6->unk_8);
        temp_s3 = (Struct_800A598C_2*)SEGMENTED_TO_VIRTUAL(phi_s6->unk_C);
        if (temp_a0_2 == 1) {
            Vec3f spAC;
            spAC.x = temp_s3[0].x;
            spAC.y = temp_s3[0].y;
            spAC.z = temp_s3[0].z;
            SkinMatrix_Vec3fMtxFMultXYZ(&D_801600C0[temp_s3[0].unk_0], &spAC, &spDC);
        } else if (arg3 == 1) {
            Struct_800A598C_2* temp_s0 = &temp_s3[phi_s6->unk_4];
            Vec3f spA0;
            spA0.x = temp_s0->x;
            spA0.y = temp_s0->y;
            spA0.z = temp_s0->z;
            SkinMatrix_Vec3fMtxFMultXYZ(&D_801600C0[temp_s0->unk_0], &spA0, &spDC);
        } else {
            Struct_800A598C_2* phi_s0 = &temp_s3[0];
            f32 phi_f26 = 0.0f;
            f32 phi_f24 = 0.0f;
            f32 phi_f22 = 0.0f;
            while (phi_s0 < &temp_s3[temp_a0_2]) {
                Vec3f sp88;
                f32 temp_f20 = phi_s0->unk_8 * 0.01f;
                sp88.x = phi_s0->x;
                sp88.y = phi_s0->y;
                sp88.z = phi_s0->z;
                SkinMatrix_Vec3fMtxFMultXYZ(&D_801600C0[phi_s0->unk_0], &sp88, &spD0);
                spD0.x *= temp_f20;
                spD0.y *= temp_f20;
                spD0.z *= temp_f20;
                phi_f26 += spD0.x;
                phi_f24 += spD0.y;
                phi_f22 += spD0.z;
                phi_s0++;
            }
            spDC.x = phi_f26;
            spDC.y = phi_f24;
            spDC.z = phi_f22;
        }
        func_800A57C0(&D_801600C0[temp_s3[phi_s6->unk_4].unk_0], temp_v0, phi_s6, spEC, &spDC);
        phi_s6++;
    }
    gSPSegment(POLY_OPA_DISP++, 0x08, avb->buf[avb->unk_0]);
    if (avb->unk_0 == 0) {
        avb->unk_0 = 1;
    } else {
        avb->unk_0 = 0;
    }

    CLOSE_DISPS(gfxCtx, "../z_skin.c", 344);
}
#else
#pragma GLOBAL_ASM("asm/non_matchings/code/z_skin/func_800A598C.s")
#endif

void func_800A5E28(GraphicsContext* gfxCtx, PSkinAwb* skin, s32 limbIndex, s32 arg3, s32 arg4) {
    SkinLimb** limbs;
    Struct_800A5E28* temp_t9;

    OPEN_DISPS(gfxCtx, "../z_skin.c", 364);
    limbs = (SkinLimb**)SEGMENTED_TO_VIRTUAL(skin->skeletonHeader->segment);
    temp_t9 = (Struct_800A5E28*)SEGMENTED_TO_VIRTUAL(((SkinLimb*)SEGMENTED_TO_VIRTUAL(limbs[limbIndex]))->segment);
    if (!(arg4 & 1)) {
        func_800A598C(gfxCtx, skin, limbIndex, arg3);
    }
    gSPDisplayList(POLY_OPA_DISP++, temp_t9->unk_8);
    CLOSE_DISPS(gfxCtx, "../z_skin.c", 377);
}

void func_800A5F60(GraphicsContext* gfxCtx, PSkinAwb* skin, s32 limbIndex, Gfx* arg3, s32 arg4) {
    Gfx* gfx = arg3;
    SkinLimb** limbs;
    s32 pad;

    OPEN_DISPS(gfxCtx, "../z_skin.c", 395);
    limbs = SEGMENTED_TO_VIRTUAL(skin->skeletonHeader->segment);
    if (arg3 == NULL) {
        gfx = (Gfx*)((SkinLimb*)SEGMENTED_TO_VIRTUAL(limbs[limbIndex]))->segment;
    }
    if (gfx != NULL) {
        Mtx* mtx = SkinMatrix_MtxFToNewMtx(gfxCtx, &D_801600C0[limbIndex]);
        if (mtx != NULL) {
            gSPMatrix(POLY_OPA_DISP++, mtx, G_MTX_PUSH | G_MTX_MUL | G_MTX_MODELVIEW);
            gSPDisplayList(POLY_OPA_DISP++, gfx);
            gSPPopMatrix(POLY_OPA_DISP++, G_MTX_MODELVIEW);
            gDPPipeSync(POLY_OPA_DISP++);
        }
    }
    CLOSE_DISPS(gfxCtx, "../z_skin.c", 433);
}

void func_800A60D8(Actor* actor, GlobalContext* globalCtx, PSkinAwb* skin, SkinCallback callback, SkinCallback2 arg4,
                   s32 arg5, s32 arg6, s32 arg7) {
    s32 i;
    s32 dataType;
    SkinLimb** limbs;
    GraphicsContext* gfxCtx = globalCtx->state.gfxCtx;
    s32 pad;

    OPEN_DISPS(gfxCtx, "../z_skin.c", 471);
    if ((arg7 & 1) == 0) {
        func_800A6AC4(skin, &D_801600C0[0], actor, arg5);
    }
    limbs = SEGMENTED_TO_VIRTUAL(skin->skeletonHeader->segment);
    if ((arg7 & 2) == 0) {
        Mtx* mtx;
        gSPMatrix(POLY_OPA_DISP++, &gMtxClear, G_MTX_NOPUSH | G_MTX_LOAD | G_MTX_MODELVIEW);
        mtx = SkinMatrix_MtxFToNewMtx(gfxCtx, &skin->mtx);
        if (mtx == NULL) {
            goto close_disps;
        }
        gSPMatrix(POLY_OPA_DISP++, mtx, G_MTX_NOPUSH | G_MTX_LOAD | G_MTX_MODELVIEW);
    }
    for (i = 0; i < skin->skeletonHeader->limbCount; i++) {
        s32 phi_a0 = 1;
        if (arg4 != NULL) {
            phi_a0 = arg4(actor, globalCtx, i, skin);
        }
        dataType = ((SkinLimb*)SEGMENTED_TO_VIRTUAL(limbs[i]))->unk_8;
        if (dataType == 4 && phi_a0 == 1) {
            func_800A5E28(gfxCtx, skin, i, arg6, arg7);
        } else if (dataType == 0xB && phi_a0 == 1) {
            func_800A5F60(gfxCtx, skin, i, NULL, arg7);
        }
    }
    if (callback != NULL) {
        callback(actor, globalCtx, skin);
    }
close_disps:
    CLOSE_DISPS(gfxCtx, "../z_skin.c", 534);
}

void func_800A6330(Actor* this, GlobalContext* globalCtx, PSkinAwb* skin, SkinCallback callback, s32 arg4) {
    func_800A60D8(this, globalCtx, skin, callback, NULL, arg4, 0, 0);
}

void func_800A6360(Actor* this, GlobalContext* globalCtx, PSkinAwb* skin, SkinCallback callback, SkinCallback2 arg4,
                   s32 arg5) {
    func_800A60D8(this, globalCtx, skin, callback, arg4, arg5, 0, 0);
}

void func_800A6394(Actor* this, GlobalContext* globalCtx, PSkinAwb* skin, SkinCallback callback, SkinCallback2 arg4,
                   s32 arg5, s32 arg6) {
    func_800A60D8(this, globalCtx, skin, callback, arg4, arg5, arg6, 0);
}

void func_800A63CC(Actor* this, GlobalContext* globalCtx, PSkinAwb* skin, SkinCallback callback, SkinCallback2 arg4,
                   s32 arg5, s32 arg6, s32 arg7) {
    func_800A60D8(this, globalCtx, skin, callback, arg4, arg5, arg6, arg7);
}

void func_800A6408(PSkinAwb* skin, s32 joint, Vec3f* arg2, Vec3f* arg3) {
    MtxF sp18;

    SkinMatrix_MtxFMtxFMult(&skin->mtx, &D_801600C0[joint], &sp18);
    SkinMatrix_Vec3fMtxFMultXYZ(&sp18, arg2, arg3);
}
