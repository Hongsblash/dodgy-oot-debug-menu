#include <ultra64.h>
#include <global.h>

u16 sReactionTextIds[60][9] = {
    { 0x0000, 0x7124, 0x7127, 0x7126, 0x7125, 0x7127, 0x7124, 0x7125, 0x7127 },
    { 0x0000, 0x7128, 0x7129, 0x7128, 0x7128, 0x7128, 0x7128, 0x712A, 0x712B },
    { 0x0000, 0x7128, 0x712B, 0x7128, 0x7128, 0x7129, 0x7128, 0x712B, 0x7128 },
    { 0x0000, 0x7128, 0x7129, 0x7128, 0x7128, 0x7128, 0x7128, 0x712A, 0x712B },
    { 0x0000, 0x7128, 0x7129, 0x712B, 0x7128, 0x7128, 0x7128, 0x7129, 0x7128 },
    { 0x0000, 0x712D, 0x712D, 0x712D, 0x712D, 0x712D, 0x712D, 0x712D, 0x712F },
    { 0x0000, 0x712C, 0x712C, 0x712C, 0x712E, 0x712C, 0x712C, 0x712F, 0x712F },
    { 0x0000, 0x712C, 0x712C, 0x712C, 0x712F, 0x712C, 0x712C, 0x712F, 0x712F },
    { 0x0000, 0x7130, 0x7132, 0x7133, 0x7130, 0x7130, 0x7131, 0x7132, 0x7131 },
    { 0x0000, 0x7134, 0x7137, 0x7135, 0x7134, 0x7136, 0x7135, 0x7134, 0x7135 },
    { 0x0000, 0x7138, 0x713A, 0x7138, 0x7139, 0x713A, 0x7138, 0x7139, 0x713B },
    { 0x0000, 0x7144, 0x7146, 0x7144, 0x7146, 0x7147, 0x7145, 0x7145, 0x7147 },
    { 0x0000, 0x7148, 0x7149, 0x7149, 0x714A, 0x714A, 0x714B, 0x7149, 0x714B },
    { 0x0000, 0x714C, 0x714D, 0x714C, 0x714C, 0x714E, 0x714C, 0x714E, 0x714F },
    { 0x0000, 0x7150, 0x7153, 0x7152, 0x7150, 0x7151, 0x7153, 0x7153, 0x7151 },
    { 0x0000, 0x7155, 0x7156, 0x7157, 0x7154, 0x7156, 0x7156, 0x7156, 0x7156 },
    { 0x0000, 0x715A, 0x7159, 0x715B, 0x715A, 0x715A, 0x7158, 0x7158, 0x715B },
    { 0x0000, 0x715E, 0x715D, 0x715D, 0x715F, 0x715E, 0x715C, 0x715C, 0x715D },
    { 0x0000, 0x7163, 0x7162, 0x7160, 0x7163, 0x7160, 0x7161, 0x7161, 0x7160 },
    { 0x0000, 0x7164, 0x7166, 0x7164, 0x7167, 0x7164, 0x7164, 0x7164, 0x7167 },
    { 0x0000, 0x716B, 0x7169, 0x7168, 0x716B, 0x716A, 0x716B, 0x716B, 0x716A },
    { 0x0000, 0x716C, 0x716D, 0x716F, 0x716C, 0x716E, 0x716E, 0x716E, 0x716F },
    { 0x0000, 0x7171, 0x7173, 0x7170, 0x7172, 0x0000, 0x0000, 0x0000, 0x0000 },
    { 0x0000, 0x7176, 0x7177, 0x7174, 0x7174, 0x7175, 0x7174, 0x7174, 0x7177 },
    { 0x0000, 0x7178, 0x7179, 0x7179, 0x717B, 0x717A, 0x717B, 0x717A, 0x717B },
    { 0x0000, 0x717D, 0x717C, 0x717C, 0x717D, 0x717F, 0x717C, 0x717E, 0x717D },
    { 0x0000, 0x7183, 0x7181, 0x7180, 0x7183, 0x7182, 0x7183, 0x7181, 0x7183 },
    { 0x0000, 0x7184, 0x7186, 0x7185, 0x7186, 0x7184, 0x7187, 0x7186, 0x7184 },
    { 0x0000, 0x71A4, 0x71A6, 0x71A5, 0x0000, 0x71A6, 0x71A6, 0x71A6, 0x71A7 },
    { 0x0000, 0x7188, 0x7188, 0x7189, 0x7188, 0x7189, 0x718B, 0x718A, 0x7189 },
    { 0x0000, 0x718C, 0x718C, 0x718D, 0x718C, 0x718E, 0x718F, 0x718D, 0x718C },
    { 0x0000, 0x7190, 0x7190, 0x7191, 0x7192, 0x7191, 0x7193, 0x7190, 0x7191 },
    { 0x0000, 0x7196, 0x7194, 0x7195, 0x7196, 0x7197, 0x7194, 0x7196, 0x7195 },
    { 0x0000, 0x7199, 0x719A, 0x7198, 0x7198, 0x719A, 0x719A, 0x719B, 0x7198 },
    { 0x0000, 0x719D, 0x719C, 0x719E, 0x719D, 0x719D, 0x719C, 0x719F, 0x719E },
    { 0x0000, 0x71A1, 0x71A0, 0x71A1, 0x71A2, 0x71A1, 0x71A2, 0x71A3, 0x71A2 },
    { 0x0000, 0x711C, 0x711E, 0x711C, 0x711F, 0x711E, 0x711C, 0x711D, 0x711F },
    { 0x0000, 0x7104, 0x7105, 0x7107, 0x7107, 0x7105, 0x7106, 0x7107, 0x7107 },
    { 0x0000, 0x7107, 0x7105, 0x7107, 0x7107, 0x7106, 0x7107, 0x7107, 0x7105 },
    { 0x0000, 0x7113, 0x7117, 0x7113, 0x7110, 0x7112, 0x7112, 0x7116, 0x7112 },
    { 0x0000, 0x7113, 0x7113, 0x7113, 0x7113, 0x7113, 0x7113, 0x7111, 0x7113 },
    { 0x0000, 0x7113, 0x7117, 0x7113, 0x7110, 0x7112, 0x7112, 0x7116, 0x7112 },
    { 0x0000, 0x7117, 0x7117, 0x7117, 0x7117, 0x7117, 0x7117, 0x7117, 0x7113 },
    { 0x0000, 0x7101, 0x7100, 0x7102, 0x7103, 0x7101, 0x7100, 0x7102, 0x7103 },
    { 0x0000, 0x7100, 0x7102, 0x7100, 0x7100, 0x7100, 0x7100, 0x7100, 0x7102 },
    { 0x0000, 0x710A, 0x7109, 0x7109, 0x710A, 0x710B, 0x7108, 0x7109, 0x710B },
    { 0x0000, 0x7117, 0x7112, 0x7113, 0x7110, 0x710C, 0x7117, 0x710E, 0x7112 },
    { 0x0000, 0x710D, 0x710F, 0x710C, 0x7112, 0x710D, 0x710C, 0x710C, 0x710F },
    { 0x0000, 0x710A, 0x7109, 0x711A, 0x710A, 0x7109, 0x7108, 0x710B, 0x7109 },
    { 0x0000, 0x710C, 0x710F, 0x7113, 0x7110, 0x710D, 0x7112, 0x7116, 0x710D },
    { 0x0000, 0x7115, 0x7114, 0x7114, 0x7115, 0x7114, 0x7114, 0x7116, 0x7117 },
    { 0x0000, 0x7113, 0x710F, 0x7113, 0x7110, 0x710C, 0x711A, 0x710D, 0x7112 },
    { 0x0000, 0x7101, 0x7102, 0x7103, 0x7101, 0x7100, 0x7100, 0x7102, 0x7100 },
    { 0x0000, 0x7112, 0x710E, 0x7112, 0x710E, 0x710D, 0x7112, 0x710E, 0x710F },
    { 0x0000, 0x7142, 0x7141, 0x7142, 0x7143, 0x7140, 0x7140, 0x7141, 0x7143 },
    { 0x0000, 0x713C, 0x713D, 0x713D, 0x713E, 0x713E, 0x713F, 0x713D, 0x713F },
    { 0x0000, 0x7101, 0x7102, 0x7103, 0x7101, 0x7100, 0x7100, 0x7102, 0x7100 },
    { 0x0000, 0x7113, 0x7117, 0x7113, 0x7110, 0x7112, 0x7112, 0x7116, 0x7112 },
    { 0x0000, 0x7104, 0x7105, 0x7107, 0x7105, 0x7105, 0x7105, 0x7107, 0x7107 },
    { 0x0000, 0x7104, 0x7105, 0x7107, 0x7105, 0x710C, 0x7105, 0x7107, 0x7107 },
};

u16 Text_GetFaceReaction(GlobalContext* globalCtx, u32 reactionSet) {
    u8 currentMask = func_8008F080(globalCtx);
    return sReactionTextIds[reactionSet][currentMask];
}
