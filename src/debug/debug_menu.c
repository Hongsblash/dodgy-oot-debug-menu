#include "global.h"
#include <ultra64.h>
#include "z64.h"

#define Align(val, align) ((((val) % (align)) != 0) ? (val) + (align) - ((val) % (align)) : (val))
#define U32_RGB(x) (u8)(x >> 24), (u8)(x >> 16), (u8)(x >> 8)

static Gfx sPolyGfxInit_HitBox[] = {
	gsSPLoadGeometryMode(G_ZBUFFER | G_SHADE | G_LIGHTING),
	gsSPTexture(0, 0, 0, G_TX_RENDERTILE, G_OFF),
	gsDPPipeSync(),
	gsDPSetCycleType(G_CYC_1CYCLE),
	gsDPSetRenderMode(
		Z_CMP | IM_RD | CVG_DST_FULL | FORCE_BL | ZMODE_XLU | GBL_c1(G_BL_CLR_IN,G_BL_A_IN, G_BL_CLR_MEM, G_BL_1MA),
		Z_CMP | IM_RD | CVG_DST_FULL | FORCE_BL | ZMODE_XLU | GBL_c2(G_BL_CLR_IN,G_BL_A_IN, G_BL_CLR_MEM, G_BL_1MA)
	),
	gsDPSetCombineLERP(PRIMITIVE, 0, SHADE, 0, 0, 0, 0, ENVIRONMENT, PRIMITIVE, 0, SHADE, 0, 0, 0, 0, ENVIRONMENT),
	gsDPSetEnvColor(255, 255, 255, 128),
	gsSPEndDisplayList(),
};

static Gfx sPolyGfxInit_Collision[] = {
	gsSPLoadGeometryMode(G_ZBUFFER | G_SHADE | G_LIGHTING),
	gsSPTexture(0, 0, 0, G_TX_RENDERTILE, G_OFF),
	gsDPPipeSync(),
	gsDPSetCycleType(G_CYC_1CYCLE),
	gsDPSetRenderMode(
		Z_CMP | IM_RD | CVG_DST_FULL | FORCE_BL | ZMODE_DEC | GBL_c1(G_BL_CLR_IN,G_BL_A_IN, G_BL_CLR_MEM, G_BL_1MA),
		Z_CMP | IM_RD | CVG_DST_FULL | FORCE_BL | ZMODE_DEC | GBL_c2(G_BL_CLR_IN,G_BL_A_IN, G_BL_CLR_MEM, G_BL_1MA)
	),
	gsDPSetCombineLERP(PRIMITIVE, 0, SHADE, 0, 0, 0, 0, ENVIRONMENT, PRIMITIVE, 0, SHADE, 0, 0, 0, 0, ENVIRONMENT),
	gsDPSetEnvColor(255, 255, 255, 128),
	gsSPEndDisplayList(),
};

#define DEBUG_GFX_BUFFER_SIZE 0x4000
static Gfx sPolyBuffer[DEBUG_GFX_BUFFER_SIZE];

void Debug_Text(PlayState* playState, u8 r, u8 g, u8 b, s32 x, s32 y, char* fmt, ...) {
    va_list args;
    GfxPrint debTex;

    OPEN_DISPS(playState->state.gfxCtx, __FILE__, __LINE__);

    Gfx* polyOpa = POLY_OPA_DISP;
    Gfx* gfx = Graph_GfxPlusOne(polyOpa);
    
    va_start(args, fmt);
    gSPDisplayList(OVERLAY_DISP++, gfx);
    GfxPrint_Init(&debTex);
    GfxPrint_Open(&debTex, gfx);
    GfxPrint_SetColor(&debTex, r, g, b, 255);
    GfxPrint_SetPos(&debTex, x, y);
    GfxPrint_VPrintf(&debTex, fmt, args);
    gfx = GfxPrint_Close(&debTex);
    GfxPrint_Destroy(&debTex);
    
    gSPEndDisplayList(gfx++);
    Graph_BranchDlist(polyOpa, gfx);
    POLY_OPA_DISP = gfx;
    
    CLOSE_DISPS(playState->state.gfxCtx, __FILE__, __LINE__);

    va_end(args);
}

void Debug_TextScaled(PlayState* playState, int scale, u8 r, u8 g, u8 b, s32 x, s32 y, char* fmt, ...) {
    va_list args;
    GfxPrint debTex;

    OPEN_DISPS(playState->state.gfxCtx, __FILE__, __LINE__);

    Gfx* polyOpa = POLY_OPA_DISP;
    Gfx* gfx = Graph_GfxPlusOne(polyOpa);
    
    va_start(args, fmt);
    gSPDisplayList(OVERLAY_DISP++, gfx);
    GfxPrint_Init(&debTex);
    GfxPrint_Open(&debTex, gfx);
    GfxPrint_SetColor(&debTex, r, g, b, 255);
    GfxPrint_SetPos(&debTex, x, y);
    GfxPrint_VPrintfScaled(&debTex, scale, fmt, args);
    gfx = GfxPrint_Close(&debTex);
    GfxPrint_Destroy(&debTex);
    
    gSPEndDisplayList(gfx++);
    Graph_BranchDlist(polyOpa, gfx);
    POLY_OPA_DISP = gfx;
    
    CLOSE_DISPS(playState->state.gfxCtx, __FILE__, __LINE__);

    va_end(args);
}

Vec3f DataMenu_TriNorm(Vec3f* v1, Vec3f* v2, Vec3f* v3) {
    Vec3f norm;
    
    norm.x = (v2->y - v1->y) * (v3->z - v1->z) - (v2->z - v1->z) * (v3->y - v1->y);
    norm.y = (v2->z - v1->z) * (v3->x - v1->x) - (v2->x - v1->x) * (v3->z - v1->z);
    norm.z = (v2->x - v1->x) * (v3->y - v1->y) - (v2->y - v1->y) * (v3->x - v1->x);
    
    f32 mag = sqrtf(SQXYZ(norm));
    
    if (mag != 0.0f) {
        norm.x *= 127.0f / mag;
        norm.y *= 127.0f / mag;
        norm.z *= 127.0f / mag;
    }
    
    return norm;
}

void DataMenu_DrawQuad3D(PlayState* playState, Gfx** gfxP, Vec3f* v1, Vec3f* v2, Vec3f* v3, Vec3f* v4) {
    Vtx* v = Graph_Alloc(playState->state.gfxCtx, 4 * sizeof(Vtx));
    Vec3f norm = DataMenu_TriNorm(v1, v2, v4);

    // Manually define vertices
    v[0].n.ob[0] = (s16) v1->x;
    v[0].n.ob[1] = (s16) v1->y;
    v[0].n.ob[2] = (s16) v1->z;
    v[0].n.flag = 0;
    v[0].n.tc[0] = 0;
    v[0].n.tc[1] = 0;
    v[0].n.n[0] = (s8) norm.x;
    v[0].n.n[1] = (s8) norm.y;
    v[0].n.n[2] = (s8) norm.z;
    v[0].n.a = 0xFF;

    v[1].n.ob[0] = (s16) v2->x;
    v[1].n.ob[1] = (s16) v2->y;
    v[1].n.ob[2] = (s16) v2->z;
    v[1].n.flag = 0;
    v[1].n.tc[0] = 0;
    v[1].n.tc[1] = 0;
    v[1].n.n[0] = (s8) norm.x;
    v[1].n.n[1] = (s8) norm.y;
    v[1].n.n[2] = (s8) norm.z;
    v[1].n.a = 0xFF;

    v[2].n.ob[0] = (s16) v3->x;
    v[2].n.ob[1] = (s16) v3->y;
    v[2].n.ob[2] = (s16) v3->z;
    v[2].n.flag = 0;
    v[2].n.tc[0] = 0;
    v[2].n.tc[1] = 0;
    v[2].n.n[0] = (s8) norm.x;
    v[2].n.n[1] = (s8) norm.y;
    v[2].n.n[2] = (s8) norm.z;
    v[2].n.a = 0xFF;

    v[3].n.ob[0] = (s16) v4->x;
    v[3].n.ob[1] = (s16) v4->y;
    v[3].n.ob[2] = (s16) v4->z;
    v[3].n.flag = 0;
    v[3].n.tc[0] = 0;
    v[3].n.tc[1] = 0;
    v[3].n.n[0] = (s8) norm.x;
    v[3].n.n[1] = (s8) norm.y;
    v[3].n.n[2] = (s8) norm.z;
    v[3].n.a = 0xFF;

    // Load vertices into the pipeline and draw the quad using two triangles
    gSPVertex((*gfxP)++, v, 4, 0);
    gSP2Triangles((*gfxP)++, 0, 1, 2, 0, 0, 2, 3, 0);
}

void DataMenu_DrawTri3D(PlayState* playState, Gfx** gfxP, Vec3f* v1, Vec3f* v2, Vec3f* v3) {
    Vtx* v = Graph_Alloc(playState->state.gfxCtx, 3 * sizeof(Vtx));
    Vec3f norm = DataMenu_TriNorm(v1, v2, v3);

    // Define vertices manually
    v[0].n.ob[0] = (s16) v1->x;
    v[0].n.ob[1] = (s16) v1->y;
    v[0].n.ob[2] = (s16) v1->z;
    v[0].n.flag = 0;
    v[0].n.tc[0] = 0;
    v[0].n.tc[1] = 0;
    v[0].n.n[0] = (s8) norm.x;
    v[0].n.n[1] = (s8) norm.y;
    v[0].n.n[2] = (s8) norm.z;
    v[0].n.a = 0xFF;

    v[1].n.ob[0] = (s16) v2->x;
    v[1].n.ob[1] = (s16) v2->y;
    v[1].n.ob[2] = (s16) v2->z;
    v[1].n.flag = 0;
    v[1].n.tc[0] = 0;
    v[1].n.tc[1] = 0;
    v[1].n.n[0] = (s8) norm.x;
    v[1].n.n[1] = (s8) norm.y;
    v[1].n.n[2] = (s8) norm.z;
    v[1].n.a = 0xFF;

    v[2].n.ob[0] = (s16) v3->x;
    v[2].n.ob[1] = (s16) v3->y;
    v[2].n.ob[2] = (s16) v3->z;
    v[2].n.flag = 0;
    v[2].n.tc[0] = 0;
    v[2].n.tc[1] = 0;
    v[2].n.n[0] = (s8) norm.x;
    v[2].n.n[1] = (s8) norm.y;
    v[2].n.n[2] = (s8) norm.z;
    v[2].n.a = 0xFF;

    // Load vertices into the pipeline and draw the triangle
    gSPVertex((*gfxP)++, v, 3, 0);
    gSP1Triangle((*gfxP)++, 0, 1, 2, 0);
}

void DataMenu_DrawCyl(PlayState* playState, Gfx** gfxP, Vec3s* pos, s16 radius, s16 height) {
    static Gfx* pCylGfx = NULL;

    if (!pCylGfx) {
        static Gfx cylGfx[5 + 12 * 2];
        static Vtx cylVtx[2 + 12 * 2];

        s32 i;
        Gfx* cylGfxP = pCylGfx = cylGfx;

        cylVtx[0].n.ob[0] = 0;
        cylVtx[0].n.ob[1] = 0;
        cylVtx[0].n.ob[2] = 0;
        cylVtx[0].n.flag = 0;
        cylVtx[0].n.tc[0] = 0;
        cylVtx[0].n.tc[1] = 0;
        cylVtx[0].n.n[0] = 0;
        cylVtx[0].n.n[1] = -127;
        cylVtx[0].n.n[2] = 0;
        cylVtx[0].n.a = 255;

        cylVtx[1].n.ob[0] = 0;
        cylVtx[1].n.ob[1] = 128;
        cylVtx[1].n.ob[2] = 0;
        cylVtx[1].n.flag = 0;
        cylVtx[1].n.tc[0] = 0;
        cylVtx[1].n.tc[1] = 0;
        cylVtx[1].n.n[0] = 0;
        cylVtx[1].n.n[1] = 127;
        cylVtx[1].n.n[2] = 0;
        cylVtx[1].n.a = 255;

        for (i = 0; i < 12; ++i) {
            s32 vtxX = Math_FFloorF(0.5f + cosf(2.f * M_PI * i / 12) * 128.f);
            s32 vtxZ = Math_FFloorF(0.5f - sinf(2.f * M_PI * i / 12) * 128.f);
            s32 normX = cosf(2.f * M_PI * i / 12) * 127.f;
            s32 normZ = -sinf(2.f * M_PI * i / 12) * 127.f;

            cylVtx[2 + i * 2 + 0].n.ob[0] = vtxX;
            cylVtx[2 + i * 2 + 0].n.ob[1] = 0;
            cylVtx[2 + i * 2 + 0].n.ob[2] = vtxZ;
            cylVtx[2 + i * 2 + 0].n.flag = 0;
            cylVtx[2 + i * 2 + 0].n.tc[0] = 0;
            cylVtx[2 + i * 2 + 0].n.tc[1] = 0;
            cylVtx[2 + i * 2 + 0].n.n[0] = normX;
            cylVtx[2 + i * 2 + 0].n.n[1] = 0;
            cylVtx[2 + i * 2 + 0].n.n[2] = normZ;
            cylVtx[2 + i * 2 + 0].n.a = 255;

            cylVtx[2 + i * 2 + 1].n.ob[0] = vtxX;
            cylVtx[2 + i * 2 + 1].n.ob[1] = 128;
            cylVtx[2 + i * 2 + 1].n.ob[2] = vtxZ;
            cylVtx[2 + i * 2 + 1].n.flag = 0;
            cylVtx[2 + i * 2 + 1].n.tc[0] = 0;
            cylVtx[2 + i * 2 + 1].n.tc[1] = 0;
            cylVtx[2 + i * 2 + 1].n.n[0] = normX;
            cylVtx[2 + i * 2 + 1].n.n[1] = 0;
            cylVtx[2 + i * 2 + 1].n.n[2] = normZ;
            cylVtx[2 + i * 2 + 1].n.a = 255;
        }

        gSPSetGeometryMode(cylGfxP++, G_CULL_BACK | G_SHADING_SMOOTH);
        gSPVertex(cylGfxP++, cylVtx, 2 + 12 * 2, 0);

        for (i = 0; i < 12; ++i) {
            s32 p = (i + 12 - 1) % 12;

            gSP2Triangles(cylGfxP++, 2 + p * 2 + 0, 2 + i * 2 + 0, 2 + i * 2 + 1, 0, 2 + p * 2 + 0, 2 + i * 2 + 1,
                          2 + p * 2 + 1, 0);
        }

        gSPClearGeometryMode(cylGfxP++, G_SHADING_SMOOTH);

        for (i = 0; i < 12; ++i) {
            s32 p = (i + 12 - 1) % 12;

            gSP2Triangles(cylGfxP++, 0, 2 + i * 2 + 0, 2 + p * 2 + 0, 0, 1, 2 + p * 2 + 1, 2 + i * 2 + 1, 0);
        }

        gSPClearGeometryMode(cylGfxP++, G_CULL_BACK);
        gSPEndDisplayList(cylGfxP++);
    }

    Matrix_Push();

    Matrix_Translate(pos->x, pos->y, pos->z, MTXMODE_NEW);
    Matrix_Scale(radius / 128.0f, height / 128.0f, radius / 128.0f, MTXMODE_APPLY);

    gSPMatrix((*gfxP)++, Matrix_NewMtx(playState->state.gfxCtx, __FILE__, __LINE__), G_MTX_MODELVIEW | G_MTX_LOAD | G_MTX_PUSH);
    gSPDisplayList((*gfxP)++, pCylGfx);
    gSPPopMatrix((*gfxP)++, G_MTX_MODELVIEW);

    Matrix_Pop();
}

void DataMenu_DrawSph(PlayState* playState, Gfx** gfxP, Vec3s* pos, s16 radius) {
    static Gfx* pSphGfx = NULL;
    static Vtx sphVtx[42];
    static Gfx sphGfx[45];

    Gfx* sphGfxP;
    s32 i;

    if (!pSphGfx) {
        Vec3f vtx[42];
        s32 r0n = 1, r0m = r0n / 5, r0i = 0 + 0;
        s32 r1n = 5, r1m = r1n / 5, r1i = r0i + r0n;
        s32 r2n = 10, r2m = r2n / 5, r2i = r1i + r1n;
        s32 r3n = 10, r3m = r3n / 5, r3i = r2i + r2n;
        s32 r4n = 10, r4m = r4n / 5, r4i = r3i + r3n;
        s32 r5n = 5, r5m = r5n / 5, r5i = r4i + r4n;
        s32 r6n = 1, r6m = r6n / 5, r6i = r5i + r5n;

        vtx[r0i + (0 * r0m + 0) % r0n].x = 0.0f;
        vtx[r0i + (0 * r0m + 0) % r0n].y = 1.0f;
        vtx[r0i + (0 * r0m + 0) % r0n].z = 0.0f;

        vtx[r6i + (0 * r6m + 0) % r6n].x = 0.0f;
        vtx[r6i + (0 * r6m + 0) % r6n].y = -1.0f;
        vtx[r6i + (0 * r6m + 0) % r6n].z = 0.0f;

        for (i = 0; i < 5; ++i) {
            static const f32 aXZ = 2.0f * M_PI / 10.0f;
            static const f32 aY = 0.463647609f; // Math_FAtanF(1.0f / 2.0f);

            vtx[r2i + (i * r2m + 0) % r2n].x = cosf(aXZ * (i * r2m + 0)) * cosf(aY * 1.0f);
            vtx[r2i + (i * r2m + 0) % r2n].y = sinf(aY * 1.0f);
            vtx[r2i + (i * r2m + 0) % r2n].z = -sinf(aXZ * (i * r2m + 0)) * cosf(aY * 1.0f);

            vtx[r4i + (i * r4m + 0) % r4n].x = cosf(aXZ * (i * r4m + 1)) * cosf(aY * -1.0f);
            vtx[r4i + (i * r4m + 0) % r4n].y = sinf(aY * -1.0f);
            vtx[r4i + (i * r4m + 0) % r4n].z = -sinf(aXZ * (i * r4m + 1)) * cosf(aY * -1.0f);
        }

        for (i = 0; i < 5; ++i) {
            Math3D_IcoSphSubdivideEdge(&vtx[r1i + (i * r1m + 0) % r1n], &vtx[r0i + (i * r0m + 0) % r0n],
                                       &vtx[r2i + (i * r2m + 0) % r2n]);
            Math3D_IcoSphSubdivideEdge(&vtx[r2i + (i * r2m + 1) % r2n], &vtx[r2i + (i * r2m + 0) % r2n],
                                       &vtx[r2i + (i * r2m + 2) % r2n]);
            Math3D_IcoSphSubdivideEdge(&vtx[r3i + (i * r3m + 0) % r3n], &vtx[r2i + (i * r2m + 0) % r2n],
                                       &vtx[r4i + (i * r4m + 0) % r4n]);
            Math3D_IcoSphSubdivideEdge(&vtx[r3i + (i * r3m + 1) % r3n], &vtx[r4i + (i * r4m + 0) % r4n],
                                       &vtx[r2i + (i * r2m + 2) % r2n]);
            Math3D_IcoSphSubdivideEdge(&vtx[r4i + (i * r4m + 1) % r4n], &vtx[r4i + (i * r4m + 0) % r4n],
                                       &vtx[r4i + (i * r4m + 2) % r4n]);
            Math3D_IcoSphSubdivideEdge(&vtx[r5i + (i * r5m + 0) % r5n], &vtx[r4i + (i * r4m + 0) % r4n],
                                       &vtx[r6i + (i * r6m + 0) % r6n]);
        }

        for (i = 0; i < 42; ++i) {
            Math3D_VtxF2L(&sphVtx[i], &vtx[i]);
        }

        sphGfxP = pSphGfx = sphGfx;

        gSPSetGeometryMode(sphGfxP++, G_CULL_BACK | G_SHADING_SMOOTH);
        gSPVertex(sphGfxP++, &sphVtx[r0i], r0n + r1n + r2n + r3n, r0i - r0i);

        r3i -= r0i;
        r2i -= r0i;
        r1i -= r0i;
        r0i -= r0i;

        for (i = 0; i < 5; ++i) {
            s32 v[24];

            v[0] = r0i + (i * r0m + 0) % r0n;
            v[1] = r1i + (i * r1m + 0) % r1n;
            v[2] = r1i + (i * r1m + 1) % r1n;
            v[3] = r1i + (i * r1m + 0) % r1n;
            v[4] = r2i + (i * r2m + 0) % r2n;
            v[5] = r2i + (i * r2m + 1) % r2n;
            v[6] = r1i + (i * r1m + 0) % r1n;
            v[7] = r2i + (i * r2m + 1) % r2n;
            v[8] = r1i + (i * r1m + 1) % r1n;
            v[9] = r1i + (i * r1m + 1) % r1n;
            v[10] = r2i + (i * r2m + 1) % r2n;
            v[11] = r2i + (i * r2m + 2) % r2n;
            v[12] = r2i + (i * r2m + 0) % r2n;
            v[13] = r3i + (i * r3m + 0) % r3n;
            v[14] = r2i + (i * r2m + 1) % r2n;
            v[15] = r2i + (i * r2m + 1) % r2n;
            v[16] = r3i + (i * r3m + 0) % r3n;
            v[17] = r3i + (i * r3m + 1) % r3n;
            v[18] = r2i + (i * r2m + 1) % r2n;
            v[19] = r3i + (i * r3m + 1) % r3n;
            v[20] = r2i + (i * r2m + 2) % r2n;
            v[21] = r2i + (i * r2m + 2) % r2n;
            v[22] = r3i + (i * r3m + 1) % r3n;
            v[23] = r3i + (i * r3m + 2) % r3n;

            gSP2Triangles(sphGfxP++, v[0], v[1], v[2], 0, v[3], v[4], v[5], 0);
            gSP2Triangles(sphGfxP++, v[6], v[7], v[8], 0, v[9], v[10], v[11], 0);
            gSP2Triangles(sphGfxP++, v[12], v[13], v[14], 0, v[15], v[16], v[17], 0);
            gSP2Triangles(sphGfxP++, v[18], v[19], v[20], 0, v[21], v[22], v[23], 0);
        }

        gSPVertex(sphGfxP++, &sphVtx[r4i], r4n + r5n + r6n, r4i - r4i);

        r6i -= r4i;
        r5i -= r4i;
        r4i -= r4i;

        for (i = 0; i < 5; ++i) {
            s32 v[24];

            v[0] = r3i + (i * r3m + 1) % r3n;
            v[1] = r4i + (i * r4m + 0) % r4n;
            v[2] = r4i + (i * r4m + 1) % r4n;
            v[3] = r3i + (i * r3m + 1) % r3n;
            v[4] = r4i + (i * r4m + 1) % r4n;
            v[5] = r3i + (i * r3m + 2) % r3n;
            v[6] = r3i + (i * r3m + 2) % r3n;
            v[7] = r4i + (i * r4m + 1) % r4n;
            v[8] = r4i + (i * r4m + 2) % r4n;
            v[9] = r3i + (i * r3m + 2) % r3n;
            v[10] = r4i + (i * r4m + 2) % r4n;
            v[11] = r3i + (i * r3m + 3) % r3n;
            v[12] = r4i + (i * r4m + 0) % r4n;
            v[13] = r5i + (i * r5m + 0) % r5n;
            v[14] = r4i + (i * r4m + 1) % r4n;
            v[15] = r4i + (i * r4m + 1) % r4n;
            v[16] = r5i + (i * r5m + 0) % r5n;
            v[17] = r5i + (i * r5m + 1) % r5n;
            v[18] = r4i + (i * r4m + 1) % r4n;
            v[19] = r5i + (i * r5m + 1) % r5n;
            v[20] = r4i + (i * r4m + 2) % r4n;
            v[21] = r5i + (i * r5m + 0) % r5n;
            v[22] = r6i + (i * r6m + 0) % r6n;
            v[23] = r5i + (i * r5m + 1) % r5n;

            gSP2Triangles(sphGfxP++, v[0], v[1], v[2], 0, v[3], v[4], v[5], 0);
            gSP2Triangles(sphGfxP++, v[6], v[7], v[8], 0, v[9], v[10], v[11], 0);
            gSP2Triangles(sphGfxP++, v[12], v[13], v[14], 0, v[15], v[16], v[17], 0);
            gSP2Triangles(sphGfxP++, v[18], v[19], v[20], 0, v[21], v[22], v[23], 0);
        }
        gSPClearGeometryMode(sphGfxP++, G_CULL_BACK | G_SHADING_SMOOTH);
        gSPEndDisplayList(sphGfxP++);
    }

    Matrix_Push();

    Matrix_Translate(pos->x, pos->y, pos->z, MTXMODE_NEW);
    Matrix_Scale(radius / 128.0f, radius / 128.0f, radius / 128.0f, MTXMODE_APPLY);

    gSPMatrix((*gfxP)++, Matrix_NewMtx(playState->state.gfxCtx, __FILE__, __LINE__), G_MTX_MODELVIEW | G_MTX_LOAD | G_MTX_PUSH);
    gSPDisplayList((*gfxP)++, pSphGfx);
    gSPPopMatrix((*gfxP)++, G_MTX_MODELVIEW);

    Matrix_Pop();
}

void DataMenu_DrawHitboxList(PlayState* playState, Gfx** gfxP, Collider** colList, s32 n, GraphicsContext* gfxCtx, u8 r, u8 g, u8 b, u8 a) {
    s32 i;

    for (i = 0; i < n; i++) {
        Collider* col = colList[i];
        
        switch (col->shape) {
            case COLSHAPE_JNTSPH: {
                ColliderJntSph* colJntSph = (ColliderJntSph*)col;
                
                s32 j;
                for (j = 0; j < colJntSph->count; j++) {
                    ColliderJntSphElement* colSphElem = &colJntSph->elements[j];
                    Sphere16 sphere = {
                        .center = colSphElem->dim.worldSphere.center,
                        .radius = colSphElem->dim.worldSphere.radius
                    };
                    DataMenu_DrawSph(
                        playState,
                        gfxP,
                        &colSphElem->dim.worldSphere.center,
                        colSphElem->dim.worldSphere.radius
                    );
                }
            }
            break;
            case COLSHAPE_CYLINDER: {
                ColliderCylinder* colCyl = (ColliderCylinder*)col;
                DataMenu_DrawCyl(
                    playState,
                    gfxP,
                    &colCyl->dim.pos,
                    colCyl->dim.radius,
                    colCyl->dim.height
                );
            }
            break;
            s32 k;
            case COLSHAPE_TRIS: {
                ColliderTris* tris = (ColliderTris*)col;
                for (k = 0; k < tris->count; k++) {
                    ColliderTrisElement* triElem = &tris->elements[i];
                    DataMenu_DrawTri3D(
                        playState,
                        gfxP,
                        &triElem->dim.vtx[0],
                        &triElem->dim.vtx[2],
                        &triElem->dim.vtx[1]
                    );
                }
            }
            break;
            case COLSHAPE_QUAD: {
                ColliderQuad* quad = (ColliderQuad*)col;
                DataMenu_DrawQuad3D(
                    playState,
                    gfxP,
                    &quad->dim.quad[0],
                    &quad->dim.quad[2],
                    &quad->dim.quad[3],
                    &quad->dim.quad[1]
                );
            }
            break;
        }
    }
}

void DataMenu_CollisionHangUp(Gfx** cur, Gfx** head, const char* message) {
    if (*cur > *head)
        Fault_AddHungupAndCrashImpl(message, "Gfx Buffer ran out of space!");
}

void DataMenu_DrawCollision(PlayState* playState, Gfx** gfxP, Gfx** gfxD, CollisionHeader* colHeader) {
    s32 prevTyid = -1;
    s16 alpha = -1;
    Camera* cam = GET_ACTIVE_CAM(playState);
    s32 i;
    s32 j;
    
    if (colHeader == NULL)
        return;
    
    for (i = 0; i < colHeader->numPolygons; i++) {
        CollisionPoly* colPoly = &colHeader->polyList[i];
        PolygonTypes* type = (PolygonTypes*)&colHeader->surfaceTypeList[colPoly->type];
        s32 tyid = -1;
        Vec3s vPos[3];
        f32 dist = MAXFLOAT;
        
        f32 currentDistance;
        for (j = 0; j < 3; j++) {
            vPos[j].x = colHeader->vtxList[colPoly->vtxData[j] & 0x1FFF].x;
            vPos[j].y = colHeader->vtxList[colPoly->vtxData[j] & 0x1FFF].y;
            vPos[j].z = colHeader->vtxList[colPoly->vtxData[j] & 0x1FFF].z;
            
            currentDistance = sqrtf((vPos[j].x - cam->eye.x) * (vPos[j].x - cam->eye.x) +
                                    (vPos[j].y - cam->eye.y) * (vPos[j].y - cam->eye.y) +
                                    (vPos[j].z - cam->eye.z) * (vPos[j].z - cam->eye.z));

            // Update the minimum distance
            if (j == 0 || currentDistance < dist) {
                dist = currentDistance;
            }
        }
        
        if (dist > 950.0f)
            continue;
        
        if (type->hookshot == 1)
            tyid = 0;
        else if (type->wallParams > SURFACE_WALL_NO_LEDGE_GRAB)
            tyid = 1;
        else if (type->floorParams == SURFACE_FLOOR_VOID)
            tyid = 2;
        else if (type->exit != 0 || type->floorParams == SURFACE_FLOOR_VOID_SMALL)
            tyid = 3;
        else if (type->behaviour != 0 || type->wallDamage)
            tyid = 4;
        else if (type->slope == 1)
            tyid = 5;
        else
            tyid = 6;
        
        if (tyid != prevTyid) {
            switch (tyid) {
                case 0:
                    gDPSetPrimColor((*gfxP)++, 0, 0, 128, 128, 255, 255);
                    break;
                case 1:
                    gDPSetPrimColor((*gfxP)++, 0, 0, 192, 0, 192, 255);
                    break;
                case 2:
                    gDPSetPrimColor((*gfxP)++, 0, 0, 255, 0, 0, 255);
                    break;
                case 3:
                    gDPSetPrimColor((*gfxP)++, 0, 0, 0, 255, 0, 255);
                    break;
                case 4:
                    gDPSetPrimColor((*gfxP)++, 0, 0, 192, 255, 192, 255);
                    break;
                case 5:
                    gDPSetPrimColor((*gfxP)++, 0, 0, 255, 255, 128, 255);
                    break;
                case 6:
                    gDPSetPrimColor((*gfxP)++, 0, 0, 255, 255, 255, 255);
                    break;
            }
            prevTyid = tyid;
        }
        
        Vtx vtx[3];
        s16 nextAlpha;
        if (dist > 100.0f) {
            f32 mul = 1.0f - (dist - 100) * 0.00111111111f;
            nextAlpha = 185 * CLAMP_MAX(mul, 1.0f);
            nextAlpha = Align(nextAlpha, 4);
            nextAlpha = CLAMP(nextAlpha, 0, 185);
        } else
            nextAlpha = 185;
        
        if (alpha != nextAlpha)
            gDPSetEnvColor((*gfxP)++, 192, 0, 192, nextAlpha);
        
        // Manually define vertices
        for (j = 0; j < 3; j++) {
            vtx[j].n.ob[0] = vPos[j].x;
            vtx[j].n.ob[1] = vPos[j].y;
            vtx[j].n.ob[2] = vPos[j].z;
            vtx[j].n.flag = 0;
            vtx[j].n.tc[0] = 0;
            vtx[j].n.tc[1] = 0;
            vtx[j].n.n[0] = colPoly->normal.x / 0x100;
            vtx[j].n.n[1] = colPoly->normal.y / 0x100;
            vtx[j].n.n[2] = colPoly->normal.z / 0x100;
            vtx[j].n.a = 0xFF;
        }
        
        Vtx* vtxPtr = Graph_Alloc(playState->state.gfxCtx, sizeof(vtx));  // Allocate memory for vertices
        
        // Copy manually defined vertices to the allocated memory
        if (vtxPtr) {
            for (j = 0; j < 3; j++) {
                vtxPtr[j].n.ob[0] = vtx[j].n.ob[0];
                vtxPtr[j].n.ob[1] = vtx[j].n.ob[1];
                vtxPtr[j].n.ob[2] = vtx[j].n.ob[2];
                vtxPtr[j].n.flag = vtx[j].n.flag;
                vtxPtr[j].n.tc[0] = vtx[j].n.tc[0];
                vtxPtr[j].n.tc[1] = vtx[j].n.tc[1];
                vtxPtr[j].n.n[0] = vtx[j].n.n[0];
                vtxPtr[j].n.n[1] = vtx[j].n.n[1];
                vtxPtr[j].n.n[2] = vtx[j].n.n[2];
                vtxPtr[j].n.a = vtx[j].n.a;
            }

            gSPVertex((*gfxP)++, vtxPtr, 3, 0);
            gSP1Triangle((*gfxP)++, 0, 1, 2, 0);
        }

        DataMenu_CollisionHangUp(gfxP, gfxD, "Scene Collision");
    }
}

void DataMenu_DrawQuad(PlayState* playState, Vec2f* pos, Vec2f* scale, Color_RGBA8* color) {
    DebugState* debugSysCtx = &sDebugMenuCtx;

    OPEN_DISPS(playState->state.gfxCtx, __FILE__, __LINE__);

    Matrix_Translate(pos->x, pos->y, 0, MTXMODE_NEW);
    Matrix_Scale(scale->x, scale->y, 1.0f, MTXMODE_APPLY);
    
    Gfx_SetupDL_39Overlay(playState->state.gfxCtx);
    Gfx_SetupDL_42Overlay(playState->state.gfxCtx);
    gSPClearGeometryMode(OVERLAY_DISP++, G_CULL_BOTH);
    
    gDPSetCombineMode(OVERLAY_DISP++, G_CC_PRIMITIVE, G_CC_PRIMITIVE);
    gDPSetPrimColor(OVERLAY_DISP++, 0, 0, color->r, color->g, color->b, color->a);
    
    gSPMatrix(OVERLAY_DISP++, Matrix_NewMtx(playState->state.gfxCtx, __FILE__, __LINE__), G_MTX_MODELVIEW | G_MTX_LOAD);
    gSPVertex(OVERLAY_DISP++, debugSysCtx->vtx, 4, 0);
    gSP1Quadrangle(OVERLAY_DISP++, 0, 2, 3, 1, 0);

    CLOSE_DISPS(playState->state.gfxCtx, __FILE__, __LINE__);
}

// PAGE FUNCS

const char* MenuStateTitles[MENU_STATE_COUNT] = {
    "Main Menu",   // Corresponds to MENU_STATE_MAIN_MENU
    "Debug Menu",  // Corresponds to MENU_STATE_DEBUG_MENU
    "UI Menu",     // Corresponds to MENU_STATE_UI_MENU
    "Weather Menu" // Corresponds to MENU_STATE_EFFECTS_MENU
};

static void DataMenu_CollisionViewUpdate(PlayState* playState) {
    DebugState* debugSysCtx = &sDebugMenuCtx;
    s32 i;
    
    if (!debugSysCtx->colViewEnabled) {
        return;
    }
    
    OPEN_DISPS(playState->state.gfxCtx, __FILE__, __LINE__);

    Gfx* dlist = sPolyBuffer;
    Gfx* pGfx = dlist;
    Gfx* dGfx = dlist + DEBUG_GFX_BUFFER_SIZE;
    
    // Setup
    gSPDisplayList(pGfx++, sPolyGfxInit_Collision);
    gSPSetGeometryMode(pGfx++, G_CULL_BACK);
    gSPMatrix(pGfx++, &gMtxClear, G_MTX_NOPUSH | G_MTX_LOAD | G_MTX_MODELVIEW);
    // Draw Static Collision
    DataMenu_DrawCollision(playState, &pGfx, &dGfx, playState->colCtx.colHeader);
    // Draw Dynapoly Collision
    for (i = 0; i < BG_ACTOR_MAX; i++) {
        if (playState->colCtx.dyna.bgActorFlags[i] & 1) {
            MtxF mtx;
            BgActor* bgActor = &playState->colCtx.dyna.bgActors[i];
            
            // Manually compute dyna SRT transformation
            Matrix_Push();
            SkinMatrix_SetTranslateRotateYXZScale(
                &mtx,
                bgActor->curTransform.scale.x,
                bgActor->curTransform.scale.y,
                bgActor->curTransform.scale.z,
                bgActor->curTransform.rot.x,
                bgActor->curTransform.rot.y,
                bgActor->curTransform.rot.z,
                bgActor->curTransform.pos.x,
                bgActor->curTransform.pos.y,
                bgActor->curTransform.pos.z
            );
            Matrix_Put(&mtx);
            gSPMatrix(
                pGfx++,
                Matrix_NewMtx(playState->state.gfxCtx, __FILE__, __LINE__),
                G_MTX_MODELVIEW | G_MTX_LOAD | G_MTX_PUSH
            );
            
            DataMenu_DrawCollision(playState, &pGfx, &dGfx, bgActor->colHeader);
            
            gSPPopMatrix(pGfx++, G_MTX_MODELVIEW);
            Matrix_Pop();
            
            DataMenu_CollisionHangUp(&pGfx, &dGfx, "BgActor Collision");
        }
    }
    // End
    gSPEndDisplayList(pGfx++);

    Debug_Text(playState, U32_RGB(0xFFFFFFFF), 1, 1, "%3.2f%c", 100 - ((f32)(dGfx - pGfx) / DEBUG_GFX_BUFFER_SIZE) * 100, '%');
    
    // Add dlist to POLY_OPA
    gSPDisplayList(POLY_OPA_DISP++, dlist);
    
    CLOSE_DISPS(playState->state.gfxCtx, __FILE__, __LINE__);
}

void DataMenu_HitboxViewUpdate(PlayState* playState, GraphicsContext* gfxCtx) {
    DebugState* debugSysCtx = &sDebugMenuCtx;
    
    if (!debugSysCtx->hitViewEnabled) {
        return;
    }
    
    OPEN_DISPS(gfxCtx, __FILE__, __LINE__);
    
    // dlist will be on POLY_OPA buffer
    Gfx* dlist = Graph_GfxPlusOne(POLY_OPA_DISP);
    Gfx* gfx = dlist;
    
    // Setup
    gSPDisplayList(gfx++, sPolyGfxInit_HitBox);
    gSPMatrix(gfx++, &gMtxClear, G_MTX_NOPUSH | G_MTX_LOAD | G_MTX_MODELVIEW);

    gDPSetPrimColor(gfx++, 0, 0, 255, 255, 255, 255);
    DataMenu_DrawHitboxList(playState, &gfx, playState->colChkCtx.colOC, playState->colChkCtx.colOCCount, gfxCtx, 255, 255, 255, 255); // OC, White
    
    gDPSetPrimColor(gfx++, 0, 0, 0, 0, 255, 255);
    DataMenu_DrawHitboxList(playState, &gfx, playState->colChkCtx.colAC, playState->colChkCtx.colACCount, gfxCtx, 0, 0, 255, 255); // AC, Blue
    
    gDPSetPrimColor(gfx++, 0, 0, 255, 0, 0, 255);
    DataMenu_DrawHitboxList(playState, &gfx, playState->colChkCtx.colAT, playState->colChkCtx.colATCount, gfxCtx, 255, 0, 0, 255); // AT, Red

    // End
    gSPEndDisplayList(gfx++);
    
    // Branch POLY_OPA past dlist
    Graph_BranchDlist(POLY_OPA_DISP, gfx);
    POLY_OPA_DISP = gfx;
    
    // Add dlist to POLY_XLU
    gSPDisplayList(POLY_XLU_DISP++, dlist);
    
    CLOSE_DISPS(gfxCtx, __FILE__, __LINE__);
}

void DataMenu_Init(DataMenu* this) {
    DebugState* debugSysCtx = &sDebugMenuCtx;
    
    static char* mainMenuItems[] = {"Debug", "UI", "Weather"};
    static char* debugMenuItems[] = {"Collision", "Hitbox", "Obj Memory", "Arena Memory", "Actor Viewer", "Gfx Buffer"};
    static char* uiMenuItems[] = {"Hide HUD", "Hide Link"};
    static char* weatherMenuItems[] = {"Start Rain", "Start Storm", "Start Snow", "lightConfig+", "lightConfig-"};

    // Assign arrays and counts to each menu state
    this->menus[MENU_STATE_MAIN_MENU] = (MenuContent){mainMenuItems, sizeof(mainMenuItems) / sizeof(char*)};
    this->menus[MENU_STATE_DEBUG_MENU] = (MenuContent){debugMenuItems, sizeof(debugMenuItems) / sizeof(char*)};
    this->menus[MENU_STATE_UI_MENU] = (MenuContent){uiMenuItems, sizeof(uiMenuItems) / sizeof(char*)};
    this->menus[MENU_STATE_EFFECTS_MENU] = (MenuContent){weatherMenuItems, sizeof(weatherMenuItems) / sizeof(char*)};

    this->menuState = MENU_STATE_MAIN_MENU;  // Start in the main menu
    this->selectedIndex = 0;  // Start at the first item

    // Initialize other properties
    debugSysCtx->isLinkHidden = false;
    debugSysCtx->isRaining = false;
    debugSysCtx->isStorming = false;
    debugSysCtx->colViewEnabled = false;
    debugSysCtx->hitViewEnabled = false;
    debugSysCtx->objMemViewEnabled = false;
    debugSysCtx->arenaMemViewEnabled = false;
    debugSysCtx->gfxBufferViewEnabled = false;
    this->dList = NULL;
    this->baseX = 0;
    this->baseY = 0;
    this->width = 0;
    this->height = 0;
    this->color.rgba = 0;
}

void DataMenu_ProcessInput(DataMenu* this, PlayState* play, GraphicsContext* gfxCtx) {
    DebugState* debugSysCtx = &sDebugMenuCtx;
    bool buttonPressed = false;
    int maxIndex = this->menus[this->menuState].itemCount - 1;

    if (CHECK_BTN_ALL(gDebug.input->press.button, BTN_DUP)) {
        if (this->selectedIndex > 0) {
            this->selectedIndex--;
            buttonPressed = true;
        }
    } else if (CHECK_BTN_ALL(gDebug.input->press.button, BTN_DDOWN)) {
        if (this->selectedIndex < maxIndex) {
            this->selectedIndex++;
            buttonPressed = true;
        }
    }

    // Process actions or changes in the menu
    if (CHECK_BTN_ALL(gDebug.input->press.button, BTN_DRIGHT)) {
        switch (this->menuState) {
            case MENU_STATE_MAIN_MENU:
                switch (this->selectedIndex) {
                    case 0:
                        this->menuState = MENU_STATE_DEBUG_MENU; // Navigate to debug menu
                        this->selectedIndex = 0;
                        break;
                    case 1:
                        this->menuState = MENU_STATE_UI_MENU;
                        this->selectedIndex = 0;
                        break;
                    case 2:
                        this->menuState = MENU_STATE_EFFECTS_MENU;
                        this->selectedIndex = 0;
                        break;
                }
                break;

            case MENU_STATE_DEBUG_MENU:
                switch (this->selectedIndex) {
                    case 0:
                        debugSysCtx->colViewEnabled = !debugSysCtx->colViewEnabled;
                        break;
                    case 1:
                        debugSysCtx->hitViewEnabled = !debugSysCtx->hitViewEnabled;
                        break;
                    case 2:
                        debugSysCtx->objMemViewEnabled = !debugSysCtx->objMemViewEnabled;
                        break;
                    case 3:
                        debugSysCtx->arenaMemViewEnabled = !debugSysCtx->arenaMemViewEnabled;
                        break;
                    case 4:
                        debugSysCtx->actorViewerEnabled = !debugSysCtx->actorViewerEnabled;
                        break;
                    case 5:
                        debugSysCtx->gfxBufferViewEnabled = !debugSysCtx->gfxBufferViewEnabled;
                        break;
                }
                break;

            case MENU_STATE_UI_MENU:
                switch (this->selectedIndex) {
                    case 0:
                        gDebug.isHudHidden = !gDebug.isHudHidden; // Toggle HUD visibility
                        break;
                    case 1:{
                        debugSysCtx->isLinkHidden = !debugSysCtx->isLinkHidden;
                        break;
                    }
                }
                break;
            case MENU_STATE_EFFECTS_MENU:
                switch (this->selectedIndex) {
                    case 0:
                        debugSysCtx->isRaining = !debugSysCtx->isRaining;

                        break;
                    case 1:{
                        debugSysCtx->isStorming = !debugSysCtx->isStorming;
                        break;
                    }
                    case 2:{
                        debugSysCtx->isSnowing = !debugSysCtx->isSnowing;
                        debugSysCtx->isActorsHidden = !debugSysCtx->isActorsHidden;
                        break;
                    }
                    case 3:{
                        play->envCtx.lightConfig ++;
                        play->envCtx.changeLightNextConfig ++;
                        break;
                    }
                    case 4:{
                        play->envCtx.lightConfig --;
                        play->envCtx.changeLightNextConfig --;
                        break;
                    }
                }
                break;
        }
        buttonPressed = true;
    }

    // Play sound if any button was pressed
    if (buttonPressed) {
        Audio_PlaySfxGeneral(NA_SE_SY_CURSOR, &gSfxDefaultPos, 4, &gSfxDefaultFreqAndVolScale, &gSfxDefaultFreqAndVolScale, &gSfxDefaultReverb);
    }
}

void DataMenu_ToggleMenu(DataMenu* this) {
    DebugState* debugSysCtx = &sDebugMenuCtx;

    // Check for BTN_DLEFT button press
    if (CHECK_BTN_ALL(gDebug.input->press.button, BTN_DLEFT)) {
        if (this->menuState > MENU_STATE_MAIN_MENU) {
            // If the menuState is above zero, pressing BTN_DLEFT will reduce the menuState, simulating a move to a previous screen
            this->menuState = MENU_STATE_MAIN_MENU;
            this->selectedIndex = 0;
        } else {
            debugSysCtx->isMenuOpen = !debugSysCtx->isMenuOpen;
        }
    }
}

// Initialize the array of actor mappings
ActorMapping actorMappings[] = {
    {ACTOR_PLAYER, "Player"},
    // {ACTOR_UNSET_1, "Unset 1"},
    // {ACTOR_EN_TEST, "En_Test"},
    // {ACTOR_UNSET_3, "Unset 3"},
    // {ACTOR_EN_GIRLA, "En_GirlA"},
    // {ACTOR_UNSET_5, "Unset 5"},
    // {ACTOR_UNSET_6, "Unset 6"},
    // {ACTOR_EN_PART, "En_Part"},
    // {ACTOR_EN_LIGHT, "En_Light"},
    // {ACTOR_EN_DOOR, "En_Door"},
    // {ACTOR_EN_BOX, "En_Box"},
    // {ACTOR_BG_DY_YOSEIZO, "Bg_Dy_Yoseizo"},
    // {ACTOR_BG_HIDAN_FIREWALL, "Bg_Hidan_Firewall"},
    // {ACTOR_EN_POH, "En_Poh"},
    // {ACTOR_EN_OKUTA, "En_Okuta"},
    // {ACTOR_BG_YDAN_SP, "Bg_Ydan_Sp"},
    // {ACTOR_EN_BOM, "En_Bom"},
    // {ACTOR_EN_WALLMAS, "En_Wallmas"},
    // {ACTOR_EN_DODONGO, "En_Dodongo"},
    // {ACTOR_EN_FIREFLY, "En_Firefly"},
    // {ACTOR_EN_HORSE, "En_Horse"},
    // {ACTOR_EN_ITEM00, "En_Item00"},
    // {ACTOR_EN_ARROW, "En_Arrow"},
    // {ACTOR_UNSET_17, "Unset 17"},
    // {ACTOR_EN_ELF, "En_Elf"},
    // {ACTOR_EN_NIW, "En_Niw"},
    // {ACTOR_UNSET_1A, "Unset 1A"},
    // {ACTOR_EN_TITE, "En_Tite"},
    // {ACTOR_EN_REEBA, "En_Reeba"},
    // {ACTOR_EN_PEEHAT, "En_Peehat"},
    // {ACTOR_EN_BUTTE, "En_Butte"},
    // {ACTOR_UNSET_1F, "Unset 1F"},
    // {ACTOR_EN_INSECT, "En_Insect"},
    // {ACTOR_EN_FISH, "En_Fish"},
    // {ACTOR_UNSET_22, "Unset 22"},
    // {ACTOR_EN_HOLL, "En_Holl"},
    // {ACTOR_EN_SCENE_CHANGE, "En_Scene_Change"},
    // {ACTOR_EN_ZF, "En_Zf"},
    // {ACTOR_EN_HATA, "En_Hata"},
    // {ACTOR_BOSS_DODONGO, "Boss_Dodongo"},
    // {ACTOR_BOSS_GOMA, "Boss_Goma"},
    // {ACTOR_EN_ZL1, "En_Zl1"},
    // {ACTOR_EN_VIEWER, "En_Viewer"},
    // {ACTOR_EN_GOMA, "En_Goma"},
    // {ACTOR_BG_PUSHBOX, "Bg_Pushbox"},
    // {ACTOR_EN_BUBBLE, "En_Bubble"},
    // {ACTOR_DOOR_SHUTTER, "Door_Shutter"},
    // {ACTOR_EN_DODOJR, "En_Dodojr"},
    // {ACTOR_EN_BDFIRE, "En_Bdfire"},
    // {ACTOR_UNSET_31, "Unset 31"},
    // {ACTOR_EN_BOOM, "En_Boom"},
    // {ACTOR_EN_TORCH2, "En_Torch2"},
    // {ACTOR_EN_BILI, "En_Bili"},
    // {ACTOR_EN_TP, "En_Tp"},
    // {ACTOR_UNSET_36, "Unset 36"},
    // {ACTOR_EN_ST, "En_St"},
    // {ACTOR_EN_BW, "En_Bw"},
    // {ACTOR_EN_A_OBJ, "En_A_Obj"},
    // {ACTOR_EN_EIYER, "En_Eiyer"},
    // {ACTOR_EN_RIVER_SOUND, "En_River_Sound"},
    // {ACTOR_EN_HORSE_NORMAL, "En_Horse_Normal"},
    // {ACTOR_EN_OSSAN, "En_Ossan"},
    // {ACTOR_BG_TREEMOUTH, "Bg_Treemouth"},
    // {ACTOR_BG_DODOAGO, "Bg_Dodoago"},
    // {ACTOR_BG_HIDAN_DALM, "Bg_Hidan_Dalm"},
    // {ACTOR_BG_HIDAN_HROCK, "Bg_Hidan_Hrock"},
    // {ACTOR_EN_HORSE_GANON, "En_Horse_Ganon"},
    // {ACTOR_BG_HIDAN_ROCK, "Bg_Hidan_Rock"},
    // {ACTOR_BG_HIDAN_RSEKIZOU, "Bg_Hidan_Rsekizou"},
    // {ACTOR_BG_HIDAN_SEKIZOU, "Bg_Hidan_Sekizou"},
    // {ACTOR_BG_HIDAN_SIMA, "Bg_Hidan_Sima"},
    // {ACTOR_BG_HIDAN_SYOKU, "Bg_Hidan_Syoku"},
    // {ACTOR_EN_XC, "En_Xc"},
    // {ACTOR_BG_HIDAN_CURTAIN, "Bg_Hidan_Curtain"},
    // {ACTOR_BG_SPOT00_HANEBASI, "Bg_Spot00_Hanebasi"},
    // {ACTOR_EN_MB, "En_Mb"},
    // {ACTOR_EN_BOMBF, "En_Bombf"},
    // {ACTOR_EN_ZL2, "En_Zl2"},
    // {ACTOR_BG_HIDAN_FSLIFT, "Bg_Hidan_Fslift"},
    // {ACTOR_EN_OE2, "En_OE2"},
    // {ACTOR_BG_YDAN_HASI, "Bg_Ydan_Hasi"},
    // {ACTOR_BG_YDAN_MARUTA, "Bg_Ydan_Maruta"},
    // {ACTOR_BOSS_GANONDROF, "Boss_Ganondrof"},
    // {ACTOR_UNSET_53, "Unset 53"},
    // {ACTOR_EN_AM, "En_Am"},
    // {ACTOR_EN_DEKUBABA, "En_Dekubaba"},
    // {ACTOR_EN_M_FIRE1, "En_M_Fire1"},
    // {ACTOR_EN_M_THUNDER, "En_M_Thunder"},
    // {ACTOR_BG_DDAN_JD, "Bg_Ddan_Jd"},
    // {ACTOR_BG_BREAKWALL, "Bg_Breakwall"},
    // {ACTOR_EN_JJ, "En_Jj"},
    // {ACTOR_EN_HORSE_ZELDA, "En_Horse_Zelda"},
    // {ACTOR_BG_DDAN_KD, "Bg_Ddan_Kd"},
    // {ACTOR_DOOR_WARP1, "Door_Warp1"},
    // {ACTOR_OBJ_SYOKUDAI, "Obj_Syokudai"},
    // {ACTOR_ITEM_B_HEART, "Item_B_Heart"},
    // {ACTOR_EN_DEKUNUTS, "En_Dekunuts"},
    // {ACTOR_BG_MENKURI_KAITEN, "Bg_Menkuri_Kaiten"},
    // {ACTOR_BG_MENKURI_EYE, "Bg_Menkuri_Eye"},
    // {ACTOR_EN_VALI, "En_Vali"},
    // {ACTOR_BG_MIZU_MOVEBG, "Bg_Mizu_Movebg"},
    // {ACTOR_BG_MIZU_WATER, "Bg_Mizu_Water"},
    // {ACTOR_ARMS_HOOK, "Arms_Hook"},
    // {ACTOR_EN_FHG, "En_fHG"},
    // {ACTOR_BG_MORI_HINERI, "Bg_Mori_Hineri"},
    // {ACTOR_EN_BB, "En_Bb"},
    // {ACTOR_BG_TOKI_HIKARI, "Bg_Toki_Hikari"},
    // {ACTOR_EN_YUKABYUN, "En_Yukabyun"},
    // {ACTOR_BG_TOKI_SWD, "Bg_Toki_Swd"},
    // {ACTOR_EN_FHG_FIRE, "En_Fhg_Fire"},
    // {ACTOR_BG_MJIN, "Bg_Mjin"},
    // {ACTOR_BG_HIDAN_KOUSI, "Bg_Hidan_Kousi"},
    // {ACTOR_DOOR_TOKI, "Door_Toki"},
    // {ACTOR_BG_HIDAN_HAMSTEP, "Bg_Hidan_Hamstep"},
    // {ACTOR_EN_BIRD, "En_Bird"},
    // {ACTOR_UNSET_73, "Unset 73"},
    // {ACTOR_UNSET_74, "Unset 74"},
    // {ACTOR_UNSET_75, "Unset 75"},
    // {ACTOR_UNSET_76, "Unset 76"},
    // {ACTOR_EN_WOOD02, "En_Wood02"},
    // {ACTOR_UNSET_78, "Unset 78"},
    // {ACTOR_UNSET_79, "Unset 79"},
    // {ACTOR_UNSET_7A, "Unset 7A"},
    // {ACTOR_UNSET_7B, "Unset 7B"},
    // {ACTOR_EN_LIGHTBOX, "En_Lightbox"},
    // {ACTOR_EN_PU_BOX, "En_Pu_box"},
    // {ACTOR_UNSET_7E, "Unset 7E"},
    // {ACTOR_UNSET_7F, "Unset 7F"},
    // {ACTOR_EN_TRAP, "En_Trap"},
    // {ACTOR_EN_AROW_TRAP, "En_Arow_Trap"},
    // {ACTOR_EN_VASE, "En_Vase"},
    // {ACTOR_UNSET_83, "Unset 83"},
    // {ACTOR_EN_TA, "En_Ta"},
    // {ACTOR_EN_TK, "En_Tk"},
    // {ACTOR_BG_MORI_BIGST, "Bg_Mori_Bigst"},
    // {ACTOR_BG_MORI_ELEVATOR, "Bg_Mori_Elevator"},
    // {ACTOR_BG_MORI_KAITENKABE, "Bg_Mori_Kaitenkabe"},
    // {ACTOR_BG_MORI_RAKKATENJO, "Bg_Mori_Rakkatenjo"},
    // {ACTOR_EN_VM, "En_Vm"},
    // {ACTOR_DEMO_EFFECT, "Demo_Effect"},
    // {ACTOR_DEMO_KANKYO, "Demo_Kankyo"},
    // {ACTOR_BG_HIDAN_FWBIG, "Bg_Hidan_Fwbig"},
    // {ACTOR_EN_FLOORMAS, "En_Floormas"},
    // {ACTOR_EN_HEISHI1, "En_Heishi1"},
    // {ACTOR_EN_RD, "En_Rd"},
    // {ACTOR_EN_PO_SISTERS, "En_Po_Sisters"},
    // {ACTOR_BG_HEAVY_BLOCK, "Bg_Heavy_Block"},
    // {ACTOR_BG_PO_EVENT, "Bg_Po_Event"},
    // {ACTOR_OBJ_MURE, "Obj_Mure"},
    // {ACTOR_EN_SW, "En_Sw"},
    // {ACTOR_BOSS_FD, "Boss_Fd"},
    // {ACTOR_OBJECT_KANKYO, "Object_Kankyo"},
    // {ACTOR_EN_DU, "En_Du"},
    // {ACTOR_EN_FD, "En_Fd"},
    // {ACTOR_EN_HORSE_LINK_CHILD, "En_Horse_Link_Child"},
    // {ACTOR_DOOR_ANA, "Door_Ana"},
    // {ACTOR_BG_SPOT02_OBJECTS, "Bg_Spot02_Objects"},
    // {ACTOR_BG_HAKA, "Bg_Haka"},
    // {ACTOR_MAGIC_WIND, "Magic_Wind"},
    // {ACTOR_MAGIC_FIRE, "Magic_Fire"},
    // {ACTOR_UNSET_A0, "Unset A0"},
    // {ACTOR_EN_RU1, "En_Ru1"},
    // {ACTOR_BOSS_FD2, "Boss_Fd2"},
    // {ACTOR_EN_FD_FIRE, "En_Fd_Fire"},
    // {ACTOR_EN_DH, "En_Dh"},
    // {ACTOR_EN_DHA, "En_Dha"},
    // {ACTOR_EN_RL, "En_Rl"},
    // {ACTOR_EN_ENCOUNT1, "En_Encount1"},
    // {ACTOR_DEMO_DU, "Demo_Du"},
    // {ACTOR_DEMO_IM, "Demo_Im"},
    // {ACTOR_DEMO_TRE_LGT, "Demo_Tre_Lgt"},
    // {ACTOR_EN_FW, "En_Fw"},
    // {ACTOR_BG_VB_SIMA, "Bg_Vb_Sima"},
    // {ACTOR_EN_VB_BALL, "En_Vb_Ball"},
    // {ACTOR_BG_HAKA_MEGANE, "Bg_Haka_Megane"},
    // {ACTOR_BG_HAKA_MEGANEBG, "Bg_Haka_MeganeBG"},
    // {ACTOR_BG_HAKA_SHIP, "Bg_Haka_Ship"},
    // {ACTOR_BG_HAKA_SGAMI, "Bg_Haka_Sgami"},
    // {ACTOR_UNSET_B2, "Unset B2"},
    // {ACTOR_EN_HEISHI2, "En_Heishi2"},
    // {ACTOR_EN_ENCOUNT2, "En_Encount2"},
    // {ACTOR_EN_FIRE_ROCK, "En_Fire_Rock"},
    // {ACTOR_EN_BROB, "En_Brob"},
    // {ACTOR_MIR_RAY, "Mir_Ray"},
    // {ACTOR_BG_SPOT09_OBJ, "Bg_Spot09_Obj"},
    // {ACTOR_BG_SPOT18_OBJ, "Bg_Spot18_Obj"},
    // {ACTOR_BOSS_VA, "Boss_Va"},
    // {ACTOR_BG_HAKA_TUBO, "Bg_Haka_Tubo"},
    // {ACTOR_BG_HAKA_TRAP, "Bg_Haka_Trap"},
    // {ACTOR_BG_HAKA_HUTA, "Bg_Haka_Huta"},
    // {ACTOR_BG_HAKA_ZOU, "Bg_Haka_Zou"},
    // {ACTOR_BG_SPOT17_FUNEN, "Bg_Spot17_Funen"},
    // {ACTOR_EN_SYATEKI_ITM, "En_Syateki_Itm"},
    // {ACTOR_EN_SYATEKI_MAN, "En_Syateki_Man"},
    // {ACTOR_EN_TANA, "En_Tana"},
    // {ACTOR_EN_NB, "En_Nb"},
    // {ACTOR_BOSS_MO, "Boss_Mo"},
    // {ACTOR_EN_SB, "En_Sb"},
    // {ACTOR_EN_BIGOKUTA, "En_Bigokuta"},
    // {ACTOR_EN_KAREBABA, "En_Karebaba"},
    // {ACTOR_BG_BDAN_OBJECTS, "Bg_Bdan_Objects"},
    // {ACTOR_DEMO_SA, "Demo_Sa"},
    // {ACTOR_DEMO_GO, "Demo_Go"},
    // {ACTOR_EN_IN, "En_In"},
    // {ACTOR_EN_TR, "En_Tr"},
    // {ACTOR_BG_SPOT16_BOMBSTONE, "Bg_Spot16_Bombstone"},
    // {ACTOR_UNSET_CE, "Unset CE"},
    // {ACTOR_BG_HIDAN_KOWARERUKABE, "Bg_Hidan_Kowarerukabe"},
    // {ACTOR_BG_BOMBWALL, "Bg_Bombwall"},
    // {ACTOR_BG_SPOT08_ICEBLOCK, "Bg_Spot08_Iceblock"},
    // {ACTOR_EN_RU2, "En_Ru2"},
    // {ACTOR_OBJ_DEKUJR, "Obj_Dekujr"},
    // {ACTOR_BG_MIZU_UZU, "Bg_Mizu_Uzu"},
    // {ACTOR_BG_SPOT06_OBJECTS, "Bg_Spot06_Objects"},
    // {ACTOR_BG_ICE_OBJECTS, "Bg_Ice_Objects"},
    // {ACTOR_BG_HAKA_WATER, "Bg_Haka_Water"},
    // {ACTOR_UNSET_D8, "Unset D8"},
    // {ACTOR_EN_MA2, "En_Ma2"},
    // {ACTOR_EN_BOM_CHU, "En_Bom_Chu"},
    // {ACTOR_EN_HORSE_GAME_CHECK, "En_Horse_Game_Check"},
    // {ACTOR_BOSS_TW, "Boss_Tw"},
    // {ACTOR_EN_RR, "En_Rr"},
    // {ACTOR_EN_BA, "En_Ba"},
    // {ACTOR_EN_BX, "En_Bx"},
    // {ACTOR_EN_ANUBICE, "En_Anubice"},
    // {ACTOR_EN_ANUBICE_FIRE, "En_Anubice_Fire"},
    // {ACTOR_BG_MORI_HASHIGO, "Bg_Mori_Hashigo"},
    // {ACTOR_BG_MORI_HASHIRA4, "Bg_Mori_Hashira4"},
    // {ACTOR_BG_MORI_IDOMIZU, "Bg_Mori_Idomizu"},
    // {ACTOR_BG_SPOT16_DOUGHNUT, "Bg_Spot16_Doughnut"},
    // {ACTOR_BG_BDAN_SWITCH, "Bg_Bdan_Switch"},
    // {ACTOR_EN_MA1, "En_Ma1"},
    // {ACTOR_BOSS_GANON, "Boss_Ganon"},
    // {ACTOR_BOSS_SST, "Boss_Sst"},
    // {ACTOR_UNSET_EA, "Unset EA"},
    // {ACTOR_UNSET_EB, "Unset EB"},
    // {ACTOR_EN_NY, "En_Ny"},
    // {ACTOR_EN_FR, "En_Fr"},
    // {ACTOR_ITEM_SHIELD, "Item_Shield"},
    // {ACTOR_BG_ICE_SHELTER, "Bg_Ice_Shelter"},
    // {ACTOR_EN_ICE_HONO, "En_Ice_Hono"},
    // {ACTOR_ITEM_OCARINA, "Item_Ocarina"},
    // {ACTOR_UNSET_F2, "Unset F2"},
    // {ACTOR_UNSET_F3, "Unset F3"},
    // {ACTOR_MAGIC_DARK, "Magic_Dark"},
    // {ACTOR_DEMO_6K, "Demo_6K"},
    // {ACTOR_EN_ANUBICE_TAG, "En_Anubice_Tag"},
    // {ACTOR_BG_HAKA_GATE, "Bg_Haka_Gate"},
    // {ACTOR_BG_SPOT15_SAKU, "Bg_Spot15_Saku"},
    // {ACTOR_BG_JYA_GOROIWA, "Bg_Jya_Goroiwa"},
    // {ACTOR_BG_JYA_ZURERUKABE, "Bg_Jya_Zurerukabe"},
    // {ACTOR_UNSET_FB, "Unset FB"},
    // {ACTOR_BG_JYA_COBRA, "Bg_Jya_Cobra"},
    // {ACTOR_BG_JYA_KANAAMI, "Bg_Jya_Kanaami"},
    // {ACTOR_FISHING, "Fishing"},
    // {ACTOR_OBJ_OSHIHIKI, "Obj_Oshihiki"},
    // {ACTOR_BG_GATE_SHUTTER, "Bg_Gate_Shutter"},
    // {ACTOR_EFF_DUST, "Eff_Dust"},
    // {ACTOR_BG_SPOT01_FUSYA, "Bg_Spot01_Fusya"},
    // {ACTOR_BG_SPOT01_IDOHASHIRA, "Bg_Spot01_Idohashira"},
    // {ACTOR_BG_SPOT01_IDOMIZU, "Bg_Spot01_Idomizu"},
    // {ACTOR_BG_PO_SYOKUDAI, "Bg_Po_Syokudai"},
    // {ACTOR_BG_GANON_OTYUKA, "Bg_Ganon_Otyuka"},
    // {ACTOR_BG_SPOT15_RRBOX, "Bg_Spot15_Rrbox"},
    // {ACTOR_BG_UMAJUMP, "Bg_Umajump"},
    // {ACTOR_UNSET_109, "Unset 109"},
    // {ACTOR_ARROW_FIRE, "Arrow_Fire"},
    // {ACTOR_ARROW_ICE, "Arrow_Ice"},
    // {ACTOR_ARROW_LIGHT, "Arrow_Light"},
    // {ACTOR_UNSET_10D, "Unset 10D"},
    // {ACTOR_UNSET_10E, "Unset 10E"},
    // {ACTOR_ITEM_ETCETERA, "Item_Etcetera"},
    // {ACTOR_OBJ_KIBAKO, "Obj_Kibako"},
    // {ACTOR_OBJ_TSUBO, "Obj_Tsubo"},
    // {ACTOR_EN_WONDER_ITEM, "En_Wonder_Item"},
    // {ACTOR_EN_IK, "En_Ik"},
    // {ACTOR_DEMO_IK, "Demo_Ik"},
    // {ACTOR_EN_SKJ, "En_Skj"},
    // {ACTOR_EN_SKJNEEDLE, "En_Skjneedle"},
    // {ACTOR_EN_G_SWITCH, "En_G_Switch"},
    // {ACTOR_DEMO_EXT, "Demo_Ext"},
    // {ACTOR_DEMO_SHD, "Demo_Shd"},
    // {ACTOR_EN_DNS, "En_Dns"},
    // {ACTOR_ELF_MSG, "Elf_Msg"},
    // {ACTOR_EN_HONOTRAP, "En_Honotrap"},
    // {ACTOR_EN_TUBO_TRAP, "En_Tubo_Trap"},
    // {ACTOR_OBJ_ICE_POLY, "Obj_Ice_Poly"},
    // {ACTOR_BG_SPOT03_TAKI, "Bg_Spot03_Taki"},
    // {ACTOR_BG_SPOT07_TAKI, "Bg_Spot07_Taki"},
    // {ACTOR_EN_FZ, "En_Fz"},
    // {ACTOR_EN_PO_RELAY, "En_Po_Relay"},
    // {ACTOR_BG_RELAY_OBJECTS, "Bg_Relay_Objects"},
    // {ACTOR_EN_DIVING_GAME, "En_Diving_Game"},
    // {ACTOR_EN_KUSA, "En_Kusa"},
    // {ACTOR_OBJ_BEAN, "Obj_Bean"},
    // {ACTOR_OBJ_BOMBIWA, "Obj_Bombiwa"},
    // {ACTOR_UNSET_128, "Unset 128"},
    // {ACTOR_UNSET_129, "Unset 129"},
    // {ACTOR_OBJ_SWITCH, "Obj_Switch"},
    // {ACTOR_OBJ_ELEVATOR, "Obj_Elevator"},
    // {ACTOR_OBJ_LIFT, "Obj_Lift"},
    // {ACTOR_OBJ_HSBLOCK, "Obj_Hsblock"},
    // {ACTOR_EN_OKARINA_TAG, "En_Okarina_Tag"},
    // {ACTOR_EN_YABUSAME_MARK, "En_Yabusame_Mark"},
    // {ACTOR_EN_GOROIWA, "En_Goroiwa"},
    // {ACTOR_EN_EX_RUPPY, "En_Ex_Ruppy"},
    // {ACTOR_EN_TORYO, "En_Toryo"},
    // {ACTOR_EN_DAIKU, "En_Daiku"},
    // {ACTOR_UNSET_134, "Unset 134"},
    // {ACTOR_EN_NWC, "En_Nwc"},
    // {ACTOR_EN_BLKOBJ, "En_Blkobj"},
    // {ACTOR_ITEM_INBOX, "Item_Inbox"},
    // {ACTOR_EN_GE1, "En_Ge1"},
    // {ACTOR_OBJ_BLOCKSTOP, "Obj_Blockstop"},
    // {ACTOR_EN_SDA, "En_Sda"},
    // {ACTOR_EN_CLEAR_TAG, "En_Clear_Tag"},
    // {ACTOR_EN_NIW_LADY, "En_Niw_Lady"},
    // {ACTOR_EN_GM, "En_Gm"},
    // {ACTOR_EN_MS, "En_Ms"},
    // {ACTOR_EN_HS, "En_Hs"},
    // {ACTOR_BG_INGATE, "Bg_Ingate"},
    // {ACTOR_EN_KANBAN, "En_Kanban"},
    // {ACTOR_EN_HEISHI3, "En_Heishi3"},
    // {ACTOR_EN_SYATEKI_NIW, "En_Syateki_Niw"},
    // {ACTOR_EN_ATTACK_NIW, "En_Attack_Niw"},
    // {ACTOR_BG_SPOT01_IDOSOKO, "Bg_Spot01_Idosoko"},
    // {ACTOR_EN_SA, "En_Sa"},
    // {ACTOR_EN_WONDER_TALK, "En_Wonder_Talk"},
    // {ACTOR_BG_GJYO_BRIDGE, "Bg_Gjyo_Bridge"},
    // {ACTOR_EN_DS, "En_Ds"},
    // {ACTOR_EN_MK, "En_Mk"},
    // {ACTOR_EN_BOM_BOWL_MAN, "En_Bom_Bowl_Man"},
    // {ACTOR_EN_BOM_BOWL_PIT, "En_Bom_Bowl_Pit"},
    // {ACTOR_EN_OWL, "En_Owl"},
    // {ACTOR_EN_ISHI, "En_Ishi"},
    // {ACTOR_OBJ_HANA, "Obj_Hana"},
    // {ACTOR_OBJ_LIGHTSWITCH, "Obj_Lightswitch"},
    // {ACTOR_OBJ_MURE2, "Obj_Mure2"},
    // {ACTOR_EN_GO, "En_Go"},
    // {ACTOR_EN_FU, "En_Fu"},
    // {ACTOR_UNSET_154, "Unset 154"},
    // {ACTOR_EN_CHANGER, "En_Changer"},
    // {ACTOR_BG_JYA_MEGAMI, "Bg_Jya_Megami"},
    // {ACTOR_BG_JYA_LIFT, "Bg_Jya_Lift"},
    // {ACTOR_BG_JYA_BIGMIRROR, "Bg_Jya_Bigmirror"},
    // {ACTOR_BG_JYA_BOMBCHUIWA, "Bg_Jya_Bombchuiwa"},
    // {ACTOR_BG_JYA_AMISHUTTER, "Bg_Jya_Amishutter"},
    // {ACTOR_BG_JYA_BOMBIWA, "Bg_Jya_Bombiwa"},
    // {ACTOR_BG_SPOT18_BASKET, "Bg_Spot18_Basket"},
    // {ACTOR_UNSET_15D, "Unset 15D"},
    // {ACTOR_EN_GANON_ORGAN, "En_Ganon_Organ"},
    // {ACTOR_EN_SIOFUKI, "En_Siofuki"},
    // {ACTOR_EN_STREAM, "En_Stream"},
    // {ACTOR_UNSET_161, "Unset 161"},
    // {ACTOR_EN_MM, "En_Mm"},
    // {ACTOR_EN_KO, "En_Ko"},
    // {ACTOR_EN_KZ, "En_Kz"},
    // {ACTOR_EN_WEATHER_TAG, "En_Weather_Tag"},
    // {ACTOR_BG_SST_FLOOR, "Bg_Sst_Floor"},
    // {ACTOR_EN_ANI, "En_Ani"},
    // {ACTOR_EN_EX_ITEM, "En_Ex_Item"},
    // {ACTOR_BG_JYA_IRONOBJ, "Bg_Jya_Ironobj"},
    // {ACTOR_EN_JS, "En_Js"},
    // {ACTOR_EN_JSJUTAN, "En_Jsjutan"},
    // {ACTOR_EN_CS, "En_Cs"},
    // {ACTOR_EN_MD, "En_Md"},
    // {ACTOR_EN_HY, "En_Hy"},
    // {ACTOR_EN_GANON_MANT, "En_Ganon_Mant"},
    // {ACTOR_EN_OKARINA_EFFECT, "En_Okarina_Effect"},
    // {ACTOR_EN_MAG, "En_Mag"},
    // {ACTOR_DOOR_GERUDO, "Door_Gerudo"},
    // {ACTOR_ELF_MSG2, "Elf_Msg2"},
    // {ACTOR_DEMO_GT, "Demo_Gt"},
    // {ACTOR_EN_PO_FIELD, "En_Po_Field"},
    // {ACTOR_EFC_ERUPC, "Efc_Erupc"},
    // {ACTOR_BG_ZG, "Bg_Zg"},
    // {ACTOR_EN_HEISHI4, "En_Heishi4"},
    // {ACTOR_EN_ZL3, "En_Zl3"},
    // {ACTOR_BOSS_GANON2, "Boss_Ganon2"},
    // {ACTOR_EN_KAKASI, "En_Kakasi"},
    // {ACTOR_EN_TAKARA_MAN, "En_Takara_Man"},
    // {ACTOR_OBJ_MAKEOSHIHIKI, "Obj_Makeoshihiki"},
    // {ACTOR_OCEFF_SPOT, "Oceff_Spot"},
    // {ACTOR_END_TITLE, "End_Title"},
    // {ACTOR_UNSET_180, "Unset 180"},
    // {ACTOR_EN_TORCH, "En_Torch"},
    // {ACTOR_DEMO_EC, "Demo_Ec"},
    // {ACTOR_SHOT_SUN, "Shot_Sun"},
    // {ACTOR_EN_DY_EXTRA, "En_Dy_Extra"},
    // {ACTOR_EN_WONDER_TALK2, "En_Wonder_Talk2"},
    // {ACTOR_EN_GE2, "En_Ge2"},
    // {ACTOR_OBJ_ROOMTIMER, "Obj_Roomtimer"},
    // {ACTOR_EN_SSH, "En_Ssh"},
    // {ACTOR_EN_STH, "En_Sth"},
    // {ACTOR_OCEFF_WIPE, "Oceff_Wipe"},
    // {ACTOR_OCEFF_STORM, "Oceff_Storm"},
    // {ACTOR_EN_WEIYER, "En_Weiyer"},
    // {ACTOR_BG_SPOT05_SOKO, "Bg_Spot05_Soko"},
    // {ACTOR_BG_JYA_1FLIFT, "Bg_Jya_1flift"},
    // {ACTOR_BG_JYA_HAHENIRON, "Bg_Jya_Haheniron"},
    // {ACTOR_BG_SPOT12_GATE, "Bg_Spot12_Gate"},
    // {ACTOR_BG_SPOT12_SAKU, "Bg_Spot12_Saku"},
    // {ACTOR_EN_HINTNUTS, "En_Hintnuts"},
    // {ACTOR_EN_NUTSBALL, "En_Nutsball"},
    // {ACTOR_BG_SPOT00_BREAK, "Bg_Spot00_Break"},
    // {ACTOR_EN_SHOPNUTS, "En_Shopnuts"},
    // {ACTOR_EN_IT, "En_It"},
    // {ACTOR_EN_GELDB, "En_GeldB"},
    // {ACTOR_OCEFF_WIPE2, "Oceff_Wipe2"},
    // {ACTOR_OCEFF_WIPE3, "Oceff_Wipe3"},
    // {ACTOR_EN_NIW_GIRL, "En_Niw_Girl"},
    // {ACTOR_EN_DOG, "En_Dog"},
    // {ACTOR_EN_SI, "En_Si"},
    // {ACTOR_BG_SPOT01_OBJECTS2, "Bg_Spot01_Objects2"},
    // {ACTOR_OBJ_COMB, "Obj_Comb"},
    // {ACTOR_BG_SPOT11_BAKUDANKABE, "Bg_Spot11_Bakudankabe"},
    // {ACTOR_OBJ_KIBAKO2, "Obj_Kibako2"},
    // {ACTOR_EN_DNT_DEMO, "En_Dnt_Demo"},
    // {ACTOR_EN_DNT_JIJI, "En_Dnt_Jiji"},
    // {ACTOR_EN_DNT_NOMAL, "En_Dnt_Nomal"},
    // {ACTOR_EN_GUEST, "En_Guest"},
    // {ACTOR_BG_BOM_GUARD, "Bg_Bom_Guard"},
    // {ACTOR_EN_HS2, "En_Hs2"},
    // {ACTOR_DEMO_KEKKAI, "Demo_Kekkai"},
    // {ACTOR_BG_SPOT08_BAKUDANKABE, "Bg_Spot08_Bakudankabe"},
    // {ACTOR_BG_SPOT17_BAKUDANKABE, "Bg_Spot17_Bakudankabe"},
    // {ACTOR_UNSET_1AA, "Unset 1AA"},
    // {ACTOR_OBJ_MURE3, "Obj_Mure3"},
    // {ACTOR_EN_TG, "En_Tg"},
    // {ACTOR_EN_MU, "En_Mu"},
    // {ACTOR_EN_GO2, "En_Go2"},
    // {ACTOR_EN_WF, "En_Wf"},
    // {ACTOR_EN_SKB, "En_Skb"},
    // {ACTOR_DEMO_GJ, "Demo_Gj"},
    // {ACTOR_DEMO_GEFF, "Demo_Geff"},
    // {ACTOR_BG_GND_FIREMEIRO, "Bg_Gnd_Firemeiro"},
    // {ACTOR_BG_GND_DARKMEIRO, "Bg_Gnd_Darkmeiro"},
    // {ACTOR_BG_GND_SOULMEIRO, "Bg_Gnd_Soulmeiro"},
    // {ACTOR_BG_GND_NISEKABE, "Bg_Gnd_Nisekabe"},
    // {ACTOR_BG_GND_ICEBLOCK, "Bg_Gnd_Iceblock"},
    // {ACTOR_EN_GB, "En_Gb"},
    // {ACTOR_EN_GS, "En_Gs"},
    // {ACTOR_BG_MIZU_BWALL, "Bg_Mizu_Bwall"},
    // {ACTOR_BG_MIZU_SHUTTER, "Bg_Mizu_Shutter"},
    // {ACTOR_EN_DAIKU_KAKARIKO, "En_Daiku_Kakariko"},
    // {ACTOR_BG_BOWL_WALL, "Bg_Bowl_Wall"},
    // {ACTOR_EN_WALL_TUBO, "En_Wall_Tubo"},
    // {ACTOR_EN_PO_DESERT, "En_Po_Desert"},
    // {ACTOR_EN_CROW, "En_Crow"},
    // {ACTOR_DOOR_KILLER, "Door_Killer"},
    // {ACTOR_BG_SPOT11_OASIS, "Bg_Spot11_Oasis"},
    // {ACTOR_BG_SPOT18_FUTA, "Bg_Spot18_Futa"},
    // {ACTOR_BG_SPOT18_SHUTTER, "Bg_Spot18_Shutter"},
    // {ACTOR_EN_MA3, "En_Ma3"},
    // {ACTOR_EN_COW, "En_Cow"},
    // {ACTOR_BG_ICE_TURARA, "Bg_Ice_Turara"},
    // {ACTOR_BG_ICE_SHUTTER, "Bg_Ice_Shutter"},
    // {ACTOR_EN_KAKASI2, "En_Kakasi2"},
    // {ACTOR_EN_KAKASI3, "En_Kakasi3"},
    // {ACTOR_OCEFF_WIPE4, "Oceff_Wipe4"},
    // {ACTOR_EN_EG, "En_Eg"},
    // {ACTOR_BG_MENKURI_NISEKABE, "Bg_Menkuri_Nisekabe"},
    // {ACTOR_EN_ZO, "En_Zo"},
    // {ACTOR_OBJ_MAKEKINSUTA, "Obj_Makekinsuta"},
    // {ACTOR_EN_GE3, "En_Ge3"},
    // {ACTOR_OBJ_TIMEBLOCK, "Obj_Timeblock"},
    // {ACTOR_OBJ_HAMISHI, "Obj_Hamishi"},
    // {ACTOR_EN_ZL4, "En_Zl4"},
    // {ACTOR_EN_MM2, "En_Mm2"},
    // {ACTOR_BG_JYA_BLOCK, "Bg_Jya_Block"},
    // {ACTOR_OBJ_WARP2BLOCK, "Obj_Warp2block"}
};

// char* findActorDescription(int id) {
//     s16 i;

//     int count = sizeof(actorMappings) / sizeof(actorMappings[0]);
//     for (i = 0; i < count; i++) {
//         if (actorMappings[i].actorId == id) {
//             return actorMappings[i].description;
//         }
//     }
//     return NULL;
// }

// void DataMenu_ActorViewer(DataMenu* this, PlayState* play) {
//     DebugState* debugSysCtx = &sDebugMenuCtx;
//     s16 actorCount = 0;
//     s16 i;

//     if (debugSysCtx->actorViewerEnabled) {
//         for (i = 0; i < ACTORCAT_MAX; i++) {
//             Actor* actor = play->actorCtx.actorLists[i].head;
//             while (actor != NULL) {
//                 actorCount++;
//                 f32 speed = actor->speed;
//                 // const char* description = findActorDescription(actor->id);

//                 // if (description) {
//                 //     Debug_TextScaled(play, 2, U32_RGB(0xFFFFFFFF), 1, 30 - actorCount,
//                 //         "ID: %s, Pos: (%.2f, %.2f, %.2f), Speed: %.2f",
//                 //         description, actor->world.pos.x, actor->world.pos.y, actor->world.pos.z, speed);
//                 // } else {
//                     Debug_TextScaled(play, actorCount, U32_RGB(0xFFFFFFFF), 1, 30 - actorCount + 2,
//                         "ID: %d, Pos: (%.2f, %.2f, %.2f), Speed: %.2f",
//                         actor->id, actor->world.pos.x, actor->world.pos.y, actor->world.pos.z, speed);
//                 // }
//             }
//         }
//     }
// }

void DataMenu_ActorViewer(DataMenu* this, PlayState* play) {
    DebugState* debugSysCtx = &sDebugMenuCtx;
    s16 actorCount;
    s16 i;

    if (debugSysCtx->actorViewerEnabled) {
        for (i = 0; i < ACTORCAT_MAX; i++) {
            Actor* actor = play->actorCtx.actorLists[i].head;
            if (actor != NULL) {
                // Ensure to check if actor is NULL before accessing its properties
                actorCount ++;
                f32 speed = actor->speed;
                s16 id = actor->id;

                Debug_TextScaled(play, 2, U32_RGB(0xFFFFFFFF), 1, 30 - actorCount - 2, 
                    "ID: %d, Pos: (%.2f, %.2f, %.2f), Speed: %.2f", 
                    id, actor->world.pos.x, actor->world.pos.y, actor->world.pos.z, speed);
            }
        }
    }
}

void DataMenu_ActorVisibility(DataMenu* this, PlayState* play) {
    DebugState* debugSysCtx = &sDebugMenuCtx;
    s16 actorCount;
    s16 i;

    for (i = 0; i < ACTORCAT_MAX; i++) {
        Actor* actor = play->actorCtx.actorLists[i].head;
        if (actor != NULL && actor->id != 0) {
            actorCount ++;
            if (debugSysCtx->isActorsHidden) {
                if (debugSysCtx->originalActorsDrawFunc[actorCount] == NULL) {
                    // Save the original draw function the first time visibility is toggled
                    debugSysCtx->originalActorsDrawFunc[actorCount] = actor->draw;
                }
                actor->draw = NULL; // Hide Link
            } else {
                if (debugSysCtx->originalActorsDrawFunc[actorCount] != NULL) {
                    // Restore the original draw function
                    actor->draw = debugSysCtx->originalActorsDrawFunc[actorCount];
                }
            }
        }
    }
}

void DataMenu_WeatherControl(PlayState* play) {
    DebugState* debugSysCtx = &sDebugMenuCtx;

    if (debugSysCtx->isRaining) {
        play->envCtx.precipitation[PRECIP_SOS_MAX] = 20;
    } else {
        play->envCtx.precipitation[PRECIP_SOS_MAX] = 0;
    }

    if (debugSysCtx->isStorming) {
        play->envCtx.stormRequest = STORM_REQUEST_START;
        play->envCtx.stormState = STORM_STATE_ON;
        play->envCtx.lightningState = LIGHTNING_ON;
    } else {
        play->envCtx.stormRequest = STORM_REQUEST_NONE;
        play->envCtx.stormState = STORM_STATE_OFF;
        play->envCtx.lightningState = LIGHTNING_OFF;
    }

    if (debugSysCtx->isSnowing) {
        // gWeatherMode == WEATHER_MODE_SNOW;
        // play->envCtx.precipitation[PRECIP_SNOW_MAX] = 64;
        // play->envCtx.precipitation[PRECIP_SNOW_CUR] = 64;
    } else {
        // play->envCtx.precipitation[PRECIP_SNOW_MAX] = 0;
        // play->envCtx.precipitation[PRECIP_SNOW_CUR] = 0;
    }
}

void DataMenu_LinkVisibility(PlayState* play) {
    DebugState* debugSysCtx = &sDebugMenuCtx;
    Player* player = GET_PLAYER(play);

    if (debugSysCtx->isLinkHidden) {
        if (debugSysCtx->originalLinkDrawFunc == NULL) {
            // Save the original draw function the first time visibility is toggled
            debugSysCtx->originalLinkDrawFunc = player->actor.draw;
        }
        player->actor.draw = NULL; // Hide Link
    } else {
        if (debugSysCtx->originalLinkDrawFunc != NULL) {
            // Restore the original draw function
            player->actor.draw = debugSysCtx->originalLinkDrawFunc;
        }
    }
}

// Main update function called each frame
void DataMenu_Update(DataMenu* this, GraphicsContext* gfxCtx, PlayState* play) {
    DebugState* debugSysCtx = &sDebugMenuCtx;
    static bool initialized = false;
    
    DataMenu_ToggleMenu(this);

    if (!initialized) {
        DataMenu_Init(this);
        initialized = true;
    }

    DataMenu_ActorVisibility(this, play);
    DataMenu_GraphicsBufferView(gfxCtx, play);
    DataMenu_ActorViewer(this, play);
    DataMenu_WeatherControl(play);
    DataMenu_LinkVisibility(play);
    DataMenu_CollisionViewUpdate(play);
    DataMenu_HitboxViewUpdate(play, gfxCtx);
    DataMenu_ZeldaArenaMemView(play);
    DataMenu_ObjectMemView(play);

    if (debugSysCtx->isMenuOpen) {
        DataMenu_DrawColorRectangle(gfxCtx); // Draw menu background
        DataMenu_ProcessInput(this, play, gfxCtx); // Process input
        DataMenu_Render(this, gfxCtx); // Render the menu
    }
}

// Render the menu
void DataMenu_Render(DataMenu* this, GraphicsContext* gfxCtx) {
    DebugState* debugSysCtx = &sDebugMenuCtx;
    s32 i;

    char** currentMenuItems = this->menus[this->menuState].items;
    int itemCount = this->menus[this->menuState].itemCount;
    char* menuTitle = MenuStateTitles[this->menuState]; 

    // Update the label for the HUD toggle based on current state
    if (this->menuState == MENU_STATE_UI_MENU) {
        currentMenuItems[0] = gDebug.isHudHidden ? "Show HUD" : "Hide HUD";
        currentMenuItems[1] = debugSysCtx->isLinkHidden ? "Show Link" : "Hide Link";
    }

    if (this->menuState == MENU_STATE_EFFECTS_MENU) {
        currentMenuItems[0] = debugSysCtx->isRaining ? "Stop Rain" : "Start Rain";
        currentMenuItems[1] = debugSysCtx->isStorming ? "Stop Storm" : "Start Storm";
        currentMenuItems[2] = debugSysCtx->isSnowing ? "Stop Snow" : "Start Snow";
    }

    OPEN_DISPS(gfxCtx, __FILE__, __LINE__);

    // Render the menu title at the top of the menu
    DataMenu_DrawText(gfxCtx, menuTitle, 3, 8, 255, 255, 255, 255);  // Draw the menu title in white

    for (i = 0; i < itemCount; i++) {
        u8 r = 175, g = 175, b = 175, a = 255;  // Default color: light grey
        if (i == this->selectedIndex) {
            r = 255; g = 255; b = 0;  // Highlight color: yellow
        }
        DataMenu_DrawText(gfxCtx, currentMenuItems[i], 3, 10 + i * 1, r, g, b, a);
    }

    CLOSE_DISPS(gfxCtx, __FILE__, __LINE__);
}

void DataMenu_DrawText(GraphicsContext* gfxCtx, const char* text, s32 x, s32 y, u8 r, u8 g, u8 b, u8 a) {
    Gfx* gfx;
    GfxPrint gfxPrint;

    OPEN_DISPS(gfxCtx, __FILE__, __LINE__);

    gfx = POLY_OPA_DISP + 1;
    gSPDisplayList(OVERLAY_DISP++, gfx);
    GfxPrint_Init(&gfxPrint);
    GfxPrint_Open(&gfxPrint, gfx);
    GfxPrint_SetColor(&gfxPrint, r, g, b, a);
    GfxPrint_SetPos(&gfxPrint, x, y);
    GfxPrint_PrintfScaled(&gfxPrint, 2, "%s", text);
    gfx = GfxPrint_Close(&gfxPrint);  
    GfxPrint_Destroy(&gfxPrint);

    gSPEndDisplayList(gfx++);
    gSPBranchList(POLY_OPA_DISP, gfx);
    POLY_OPA_DISP = gfx;

    CLOSE_DISPS(gfxCtx, __FILE__, __LINE__);
}

#define BLANK 0, 0, 0, ENVIRONMENT, 0, 0, 0, ENVIRONMENT

void DataMenu_DrawColorRectangle(GraphicsContext* gfxCtx) {
    Color_RGBA8 rgba = {10, 10, 10, 200};  // RGBA color: Grey with some transparency
    Gfx* gfx;
    s32 x1 = 20;  // Left x-coordinate
    s32 y1 = 60;  // Top y-coordinate
    s32 x2 = 90; // Right x-coordinate
    s32 y2 = 130;  // Bottom y-coordinate

    // Correct and clamp coordinates
    if (x2 < x1) { s32 temp = x2; x2 = x1; x1 = temp; }
    if (y2 < y1) { s32 temp = y2; y2 = y1; y1 = temp; }
    x1 = MAX(x1, 0); y1 = MAX(y1, 0);
    x2 = MIN(x2, SCREEN_WIDTH); y2 = MIN(y2, SCREEN_HEIGHT);

    OPEN_DISPS(gfxCtx, __BASE_FILE__, __LINE__);

    // Set render and combine modes for the rectangle
    gDPSetRenderMode(OVERLAY_DISP++, G_RM_XLU_SURF, G_RM_XLU_SURF2);
    gDPSetCombineMode(OVERLAY_DISP++, BLANK, BLANK);
    gSPClearGeometryMode(OVERLAY_DISP++, G_ZBUFFER);
    gDPSetDepthSource(OVERLAY_DISP++, G_ZS_PIXEL);

    // Set fill color and draw the rectangle
    gDPSetPrimColor(OVERLAY_DISP++, 0, 0, rgba.r, rgba.g, rgba.b, rgba.a);
    gDPSetEnvColor(OVERLAY_DISP++, rgba.r, rgba.g, rgba.b, rgba.a);
    gDPFillRectangle(OVERLAY_DISP++, x1, y1, x2, y2);

    gDPPipeSync(OVERLAY_DISP++);

    // Restore depth buffer state
    gSPSetGeometryMode(OVERLAY_DISP++, G_ZBUFFER);

    CLOSE_DISPS(gfxCtx, __BASE_FILE__, __LINE__);
}

void DataMenu_Destroy(DataMenu* this) {

}

u32 CalculateGfxBufferUsage(TwoHeadGfxArena* arena) {
    if (arena->p && arena->start) {
        return (u32)((char*)arena->p - (char*)arena->start);
    }
    return 0;
}

void InitializeQuadVertices(Vtx* vtx, s16 x, s16 y, s16 width, s16 height, s16 u, s16 v) {
    vtx[0].n.ob[0] = x;
    vtx[0].n.ob[1] = y + height;
    vtx[0].n.ob[2] = 0;
    vtx[0].n.flag = 0;
    vtx[0].n.tc[0] = 0;
    vtx[0].n.tc[1] = 0;
    vtx[0].n.n[0] = 0;
    vtx[0].n.n[1] = 0;
    vtx[0].n.n[2] = 127;
    vtx[0].n.a = 255;

    vtx[1].n.ob[0] = x + width;
    vtx[1].n.ob[1] = y + height;
    vtx[1].n.ob[2] = 0;
    vtx[1].n.flag = 0;
    vtx[1].n.tc[0] = u;
    vtx[1].n.tc[1] = 0;
    vtx[1].n.n[0] = 0;
    vtx[1].n.n[1] = 0;
    vtx[1].n.n[2] = 127;
    vtx[1].n.a = 255;

    vtx[2].n.ob[0] = x;
    vtx[2].n.ob[1] = y;
    vtx[2].n.ob[2] = 0;
    vtx[2].n.flag = 0;
    vtx[2].n.tc[0] = 0;
    vtx[2].n.tc[1] = v;
    vtx[2].n.n[0] = 0;
    vtx[2].n.n[1] = 0;
    vtx[2].n.n[2] = 127;
    vtx[2].n.a = 255;

    vtx[3].n.ob[0] = x + width;
    vtx[3].n.ob[1] = y;
    vtx[3].n.ob[2] = 0;
    vtx[3].n.flag = 0;
    vtx[3].n.tc[0] = u;
    vtx[3].n.tc[1] = v;
    vtx[3].n.n[0] = 0;
    vtx[3].n.n[1] = 0;
    vtx[3].n.n[2] = 127;
    vtx[3].n.a = 255;
}

void DataMenu_GraphicsBufferView(GraphicsContext* gfxCtx, PlayState* playState) {
    DebugState* debugSysCtx = &sDebugMenuCtx;
    s16 yOffset = 20;
    s16 i;
    Color_RGBA8 memColor[] = {
        { 240, 104,  77,    0xD8 },
        { 240, 169,  77,    0xD8 },
        { 191, 240,  77,    0xD8 },
        { 77,  240,  153,   0xD8 },
        { 77,  205,  240,   0xD8 },
        { 77,  115,  240,   0xD8 },
        { 137, 77,   240,   0xD8 },
        { 240, 77,   175,   0xD8 },
    };

    if (!debugSysCtx->gfxBufferViewEnabled) {
        return;
    }

    // Title
    Debug_TextScaled(playState, 2, U32_RGB(0x146EFFFF), 1, yOffset + 2, "Graphics Buffer Usage", 0);

    // Manually create vertices for a quad
    Vtx* vtx = Graph_Alloc(gfxCtx, 4 * sizeof(Vtx));
    s16 x = 0, y = -5, width = 1, height = 10, u = 16, v = 16;
    InitializeQuadVertices(vtx, x, y, width, height, u, v);
    debugSysCtx->vtx = vtx;
    
    // Positioning and scaling of the 2 quads that make up the background of the bar
    Vec2f pos = { -150.0f, 55.0f - (yOffset * 8)};
    Vec2f scale = { 300.0f, 1.10f };
    Color_RGBA8 color = { 0x6E, 0x14, 0xFF, 0xFF };
    
    scale.x += 2.0f;
    pos.x -= 1.0f;
    scale.y *= 1.2f;
    color = (Color_RGBA8) { 0x6E, 0x14, 0xFF, 0xFF };

    scale.x = (scale.x / 2.0f);
    scale.y = (scale.y / 2.0f);

    for (i = 0; i < 3; i++) {
        DataMenu_DrawQuad(playState, &pos, &scale, &color);
        pos.y = pos.y + 16;
    }

    pos.y = 55.0f - (yOffset * 8);
    pos.x = -150.0f;
    scale = (Vec2f) { 300.0f, 1.10f };
    color = (Color_RGBA8) { 0x30, 0x30, 0x30, 0xFF };

    scale.x = (scale.x / 2.0f) - 1.0f;
    scale.y = (scale.y / 2.0f);

    for (i = 0; i < 3; i++) {
        DataMenu_DrawQuad(playState, &pos, &scale, &color);
        pos.y = pos.y + 16;
    }

    pos.y = pos.y - 16;

    // Calculate and fill the used portions of memory onto the background of the bar, 
    // to give off the appearance of a gauge
    for (i = 0; i < 3; i++) {
        float usageKB = 0;
        float totalKB = 0;
        char* bufferName = "";

        switch (i) {
            case 0: 
                usageKB = CalculateGfxBufferUsage(&gfxCtx->polyOpa) / 1024.0f;
                totalKB = sizeof(((GfxPool*)0)->polyOpaBuffer) / 1024.0f;
                bufferName = "Poly Opa";
                break;
            case 1: 
                usageKB = CalculateGfxBufferUsage(&gfxCtx->polyXlu) / 1024.0f;
                totalKB = sizeof(((GfxPool*)0)->polyXluBuffer) / 1024.0f;
                bufferName = "Poly Xlu";
                break;
            case 2: 
                usageKB = CalculateGfxBufferUsage(&gfxCtx->overlay) / 1024.0f;
                totalKB = sizeof(((GfxPool*)0)->overlayBuffer) / 1024.0f;
                bufferName = "Overlay";
                break;
        }

        // Only proceed if there is usage to display
        if (usageKB > 0) {
            // Display used and total for each buffer in KB
            Debug_TextScaled(playState, 2, U32_RGB(0xFFFFFF00), 1, (yOffset + 3) + 2 * i,
                "%s: %.2f KB / %.2f KB used", bufferName, usageKB, totalKB);

            float filledRatio = usageKB / totalKB;
            pos.x = -150.0f; // Reset starting position
            scale.x = filledRatio * (300.0f / 2); // Calculate width based on usage

            color = memColor[i % 8];
            color.a = 0xE0; // Set alpha

            DataMenu_DrawQuad(playState, &pos, &scale, &color);

            pos.y = pos.y - 16;
        }
    }
}


void DataMenu_ObjectMemView(PlayState* playState) {
    DebugState* debugSysCtx = &sDebugMenuCtx;
    static u32 selectedMem = 0;
    s32 i;
    u32 selectedColor[] = {
        0x404040FF,
        0xFFFFFF00,
    };
    Color_RGBA8 memColor[] = {
        { 240, 104,  77,    0xF8 - 0x20    },
        { 240, 169,  77,    0xF8 - 0x20    },
        { 191, 240,  77,    0xF8 - 0x20    },
        { 77,  240,  153,   0xF8 - 0x20    },
        { 77,  205,  240,   0xF8 - 0x20    },
        { 77,  115,  240,   0xF8 - 0x20    },
        { 137, 77,   240,   0xF8 - 0x20    },
        { 240, 77,   175,   0xF8 - 0x20    },
    };

    s16 yOffset = 20;

    if (!debugSysCtx->objMemViewEnabled) {
        return;
    }

    Debug_Text(
        playState,
        U32_RGB(0x146EFFFF),
        1,
        5 + yOffset,
        "Object Memory View",
        0
    );
    
    // Manually create vertices for a quad
    Vtx* vtx = Graph_Alloc(playState->state.gfxCtx, 4 * sizeof(Vtx));
    s16 x = 0, y = -5, width = 1, height = 10, u = 16, v = 16;

    InitializeQuadVertices(vtx, x, y, width, height, u, v);

    debugSysCtx->vtx = vtx;
    
    Vec2f pos = { -150.0f, 55.0f - (yOffset * 8)};
    Vec2f scale = { 300.0f, 1.10f };
    Color_RGBA8 color;
    
    scale.x += 2.0f;
    pos.x -= 1.0f;
    scale.y *= 1.2f;
    color = (Color_RGBA8) { 0x6E, 0x14, 0xFF, 0xFF };
    DataMenu_DrawQuad(playState, &pos, &scale, &color);
    
    pos.x = -150.0f;
    scale = (Vec2f) { 300.0f, 1.10f };
    color = (Color_RGBA8) { 0x30, 0x30, 0x30, 0xFF };
    DataMenu_DrawQuad(playState, &pos, &scale, &color);
    
    ObjectContext* objCtx = &playState->objectCtx;
    u32 objMemSize = (u32)objCtx->spaceEnd - (u32)objCtx->spaceStart;
    
    Debug_Text(
        playState,
        U32_RGB(selectedColor[selectedMem == 0]),
        1,
        6 + yOffset,
        "Size: %0.2f / %0.2fMB",
        (f32)((u32)objCtx->slots[objCtx->numEntries].segment - (u32)objCtx->spaceStart) / 0x100000,
        (f32)objMemSize / 0x100000
    );
    
    if (objMemSize == 0 || objCtx->numEntries == 0)
        return;
    
    for (i = 0; i <= objCtx->numEntries; i++) {
        if (objCtx->slots[i].id <= 0)
            continue;
        
        u32 objSize = gObjectTable[objCtx->slots[i].id].vromEnd - gObjectTable[objCtx->slots[i].id].vromStart;
        
        pos.x = (f32)((u32)objCtx->slots[i].segment - (u32)objCtx->spaceStart) / objMemSize;
        pos.x *= 300.0f;
        scale.x = (f32)objSize / objMemSize;
        
        color = memColor[i % 8];
        color.a = 0xE0;
        
        pos.x -= 150.0f;
        scale.x *= 300.0f;
        DataMenu_DrawQuad(playState, &pos, &scale, &color);
    }
}

void DataMenu_ZeldaArenaMemView(PlayState* playState) {
    DebugState* debugSysCtx = &sDebugMenuCtx;
    static u32 selectedMem = 0;
    u32 selectedColor[] = {
        0x404040FF,
        0xFFFFFF00,
    };
    Color_RGBA8 memColor[] = {
        { 240, 104,  77,    0xF8 - 0x20    },
        { 240, 169,  77,    0xF8 - 0x20    },
        { 191, 240,  77,    0xF8 - 0x20    },
        { 77,  240,  153,   0xF8 - 0x20    },
        { 77,  205,  240,   0xF8 - 0x20    },
        { 77,  115,  240,   0xF8 - 0x20    },
        { 137, 77,   240,   0xF8 - 0x20    },
        { 240, 77,   175,   0xF8 - 0x20    },
    };

    s16 yOffset = 16;

    if (!debugSysCtx->arenaMemViewEnabled) {
        return;
    }

    Debug_Text(
        playState,
        U32_RGB(0x146EFFFF),
        1,
        5 + yOffset,
        "ZeldaArena Memory View",
        0
    );
    
    // Manually create vertices for a quad
    Vtx* vtx = Graph_Alloc(playState->state.gfxCtx, 4 * sizeof(Vtx));
    s16 x = 0, y = -5, width = 1, height = 10, u = 16, v = 16;

    InitializeQuadVertices(vtx, x, y, width, height, u, v);

    debugSysCtx->vtx = vtx;    
    Vec2f pos = { -150.0f, 55.0f - (yOffset * 8)};
    Vec2f scale = { 300.0f, 1.10f };
    Color_RGBA8 color;
    
    scale.x += 2.0f;
    pos.x -= 1.0f;
    scale.y *= 1.2f;
    color = (Color_RGBA8) { 0x6E, 0x14, 0xFF, 0xFF };
    DataMenu_DrawQuad(playState, &pos, &scale, &color);
    
    pos.x = -150.0f;
    scale = (Vec2f) { 300.0f, 1.10f };
    color = (Color_RGBA8) { 0x30, 0x30, 0x30, 0xFF };
    DataMenu_DrawQuad(playState, &pos, &scale, &color);
    
    u32 outMaxFree;
    u32 outFree;
    u32 outAlloc;
    
    ZeldaArena_GetSizes(&outMaxFree, &outFree, &outAlloc);
    
    Debug_Text(
        playState,
        U32_RGB(selectedColor[selectedMem == 0]),
        1,
        6 + yOffset,
        "Size: %0.2f / %0.2fMB",
        (f32)((outAlloc + outFree) - outFree) / 0x100000,
        (f32)(outAlloc + outFree) / 0x100000
    );
    
    ArenaNode* arena = sZeldaArena.head;
    u32 i = 0;
    
    if (!arena)
        return;
    
    while (arena) {
        f32 arenaSize = arena->size;
        f32 arenaStart = (u32)arena - (u32)sZeldaArena.head;
        
        pos.x = arenaStart / (outAlloc + outFree);
        pos.x *= 300.0f;
        scale.x = arenaSize / (outAlloc + outFree);
        
        color = memColor[i % 8];
        if (arena->isFree) {
            color = (Color_RGBA8) { 0x40, 0x40, 0x40, 0xF8 - 0x20 };
        }
        color.a = 0xE0;
        
        pos.x -= 150.0f;
        scale.x *= 300.0f;
        DataMenu_DrawQuad(playState, &pos, &scale, &color);
        i++;
        
        arena = arena->next;
    }
}

