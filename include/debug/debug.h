#ifndef DEBUG_H
#define DEBUG_H

#include "ultra64.h"
#include "macros.h"
#include "padmgr.h"
#include "color.h"

#define COLOR_WHITE (0xFFFFFF)
#define COLOR_RED (0xFF0000)
#define COLOR_GREEN (0x00FF00)
#define COLOR_BLUE (0x0000FF)
#define COLOR_BLUE2 (0x0080FF)
#define COLOR_BLUE3 (0x00BFFF)

typedef void (* PageFunc)(PlayState*);

typedef enum {
    DEBUGSYS_PAGE_MAIN,
    DEBUGSYS_PAGE_MINIMAP
} DebugPage;

typedef struct {
    char** items;  // Pointer to an array of strings
    int itemCount;  // Number of items
} MenuContent;

typedef enum {
    MENU_STATE_MAIN_MENU,
    MENU_STATE_DEBUG_MENU,
    MENU_STATE_UI_MENU,
    MENU_STATE_EFFECTS_MENU,
    MENU_STATE_COUNT  // This keeps track of the number of states
} MenuState;

typedef struct {
    struct {
        u8 on : 1;
    } state;
    DebugPage page;
    DebugPage setPage;
    Vtx*      vtx;
    bool      isMenuOpen;
    ActorFunc originalLinkDrawFunc;
    ActorFunc originalActorsDrawFunc[16];
    bool      isActorsHidden;
    bool      isLinkHidden;
    bool      isRaining;
    bool      isStorming;
    bool      isSnowing;
    bool      actorViewerEnabled;
    bool      hitViewEnabled;
    bool      colViewEnabled;
    bool      objMemViewEnabled;
    bool      gfxBufferViewEnabled;
    bool      arenaMemViewEnabled;
    bool      cineModeEnabled;
} DebugState;

typedef struct {
    PageFunc func;
    char*    name;
    u8       playerFreeze : 2;
    u8       toggle       : 2;
    u8       state        : 2;
} DebugPageInfo;

extern Arena sZeldaArena;
static DebugState sDebugMenuCtx;

#define MENU_ITEM_COUNT 15

typedef struct {
    GameState state;
    Gfx* dList;  // Internal display list storage
    int dListIndex; // Index to track display list commands
    u32 baseX;
    u32 baseY;
    u32 width;
    u32 height;
    s32 selectedIndex;      // Index of the currently selected menu item
    Color_RGBA8_u32 color;
    MenuState menuState;
    MenuContent menus[4];  // Supports three menu states
    bool hitboxViewEnabled;
} DataMenu;

struct PlayState;

typedef struct Debug {
    PlayState* play;
    Input* input;
    GraphicsContext* gfxCtx;
    bool      isHudHidden;
    // PrintUtils printer;
// #ifdef ENABLE_INV_EDITOR
//     InventoryDebug invDebug;
// #endif
} Debug;

typedef struct {
    const char* actorId;  // This should be a string that corresponds to the ID
    const char* description;  // Description or name of the actor
} ActorMapping;

typedef struct {
    // HIGH
    struct {
        u32 blockEpona   : 1;
        u32 lowerSurface : 1;
        u32 floorParams  : 4;
        u32 wallParams   : 5;
        u32 unk_001C0000 : 3;
        u32 behaviour    : 5;
        u32 exit         : 5;
        u32 camDataIndex : 8;
    };
    // LOW
    struct {
        u32 pad           : 4;
        u32 wallDamage    : 1;
        u32 conveyorDir   : 6;
        u32 conveyorSpeed : 3;
        u32 hookshot      : 1;
        u32 echo          : 6;
        u32 lightParams   : 5;
        u32 slope         : 2;
        u32 sfx           : 4;
    };
} PolygonTypes;

typedef enum {
    SURFACE_FLOOR_VOID_SMALL        = 0x5,
    SURFACE_FLOOR_HANG_LEDGE        = 0x6,
    SURFACE_FLOOR_STOP_AIR_MOMENTUM = 0x8,
    SURFACE_FLOOR_NO_LEDGE_JUMP     = 0x9,
    SURFACE_FLOOR_DIVE              = 0xB,
    SURFACE_FLOOR_VOID              = 0xC,
} SurfaceFloorParams;

typedef enum {
    SURFACE_WALL_NO_LEDGE_GRAB = 0x1,
    SURFACE_WALL_LADDER        = 0x2,
    SURFACE_WALL_LADDER_TOP    = 0x3,
    SURFACE_WALL_VINE          = 0x4,
    SURFACE_WALL_CRAWL_A       = 0x5,
    SURFACE_WALL_CRAWL_B       = 0x6,
    SURFACE_WALL_PUSH          = 0x7,
} SurfaceWallParams;

typedef enum {
    SURFACE_BEHAVIOUR_SPECIAL_UNK_1 = 0x1,
    SURFACE_BEHAVIOUR_HURT_SPIKES,
    SURFACE_BEHAVIOUR_HURT_LAVA,
    SURFACE_BEHAVIOUR_SAND,
    SURFACE_BEHAVIOUR_SLIPPERY,
    SURFACE_BEHAVIOUR_NO_FALL_DAMAGE,
    SURFACE_BEHAVIOUR_QUICKSAND,
    SURFACE_BEHAVIOUR_JABU_WALL,
    SURFACE_BEHAVIOUR_VOID_ON_CONTACT,
    SURFACE_BEHAVIOUR_UNK_A,
    SURFACE_BEHAVIOUR_LINK_LOOK_UP,
    SURFACE_BEHAVIOUR_QUICKSAND_EPONA
} SurfaceSpecialParams;

#endif // DEBUG_H
