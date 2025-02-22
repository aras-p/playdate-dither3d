// SPDX-License-Identifier: Unlicense

#include "mathlib.h"

float3 v3_tri_normal(const float3* p1, const float3* p2, const float3* p3)
{
	float3 v = v3_cross(
		(float3) { p2->x - p1->x, p2->y - p1->y, p2->z - p1->z },
		(float3) { p3->x - p1->x, p3->y - p1->y, p3->z - p1->z });
	return v3_normalize(v);
}

xform xform_identity = { .m = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, .x = 0, .y = 0, .z = 0 };

xform xform_make(float m11, float m12, float m13, float m21, float m22, float m23, float m31, float m32, float m33)
{
	return (xform) { .m = { {m11, m12, m13}, {m21, m22, m23}, {m31, m32, m33} }, .x = 0, .y = 0, .z = 0 };
}

xform xform_axes(float3 ax, float3 ay, float3 az, float3 pos)
{
	return (xform) { .m = { {ax.x, ax.y, ax.z}, {ay.x, ay.y, ay.z}, {az.x, az.y, az.z} }, .x = pos.x, .y = pos.y, .z = pos.z };
}

xform xform_make_translate(float dx, float dy, float dz)
{
	return (xform) { .m = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} }, .x = dx, .y = dy, .z = dz };
}

xform xform_make_axis_angle(float angle, float3 axis)
{
	xform p;
	float c = cosf(angle);
	float s = sinf(angle);
	float x = axis.x;
	float y = axis.y;
	float z = axis.z;
	float d = sqrtf(x * x + y * y + z * z);
	x /= d;
	y /= d;
	z /= d;
	p.m[0][0] = c + x * x * (1 - c);
	p.m[0][1] = x * y * (1 - c) - z * s;
	p.m[0][2] = x * z * (1 - c) + y * s;
	p.m[1][0] = y * x * (1 - c) + z * s;
	p.m[1][1] = c + y * y * (1 - c);
	p.m[1][2] = y * z * (1 - c) - x * s;
	p.m[2][0] = z * x * (1 - c) - y * s;
	p.m[2][1] = z * y * (1 - c) + x * s;
	p.m[2][2] = c + z * z * (1 - c);
	p.x = p.y = p.z = 0;
	return p;
}

xform xform_multiply(const xform* l, const xform* r)
{
	xform m;
	m.m[0][0] = l->m[0][0] * r->m[0][0] + l->m[1][0] * r->m[0][1] + l->m[2][0] * r->m[0][2];
	m.m[1][0] = l->m[0][0] * r->m[1][0] + l->m[1][0] * r->m[1][1] + l->m[2][0] * r->m[1][2];
	m.m[2][0] = l->m[0][0] * r->m[2][0] + l->m[1][0] * r->m[2][1] + l->m[2][0] * r->m[2][2];
	m.m[0][1] = l->m[0][1] * r->m[0][0] + l->m[1][1] * r->m[0][1] + l->m[2][1] * r->m[0][2];
	m.m[1][1] = l->m[0][1] * r->m[1][0] + l->m[1][1] * r->m[1][1] + l->m[2][1] * r->m[1][2];
	m.m[2][1] = l->m[0][1] * r->m[2][0] + l->m[1][1] * r->m[2][1] + l->m[2][1] * r->m[2][2];
	m.m[0][2] = l->m[0][2] * r->m[0][0] + l->m[1][2] * r->m[0][1] + l->m[2][2] * r->m[0][2];
	m.m[1][2] = l->m[0][2] * r->m[1][0] + l->m[1][2] * r->m[1][1] + l->m[2][2] * r->m[1][2];
	m.m[2][2] = l->m[0][2] * r->m[2][0] + l->m[1][2] * r->m[2][1] + l->m[2][2] * r->m[2][2];
	m.x = l->x * r->m[0][0] + l->y * r->m[0][1] + l->z * r->m[0][2] + r->x;
	m.y = l->x * r->m[1][0] + l->y * r->m[1][1] + l->z * r->m[1][2] + r->y;
	m.z = l->x * r->m[2][0] + l->y * r->m[2][1] + l->z * r->m[2][2] + r->z;
	return m;
}

float3 xform_transform_pt(const xform* m, float3 p)
{
	float x = p.x * m->m[0][0] + p.y * m->m[0][1] + p.z * m->m[0][2] + m->x;
	float y = p.x * m->m[1][0] + p.y * m->m[1][1] + p.z * m->m[1][2] + m->y;
	float z = p.x * m->m[2][0] + p.y * m->m[2][1] + p.z * m->m[2][2] + m->z;

	return (float3) { x, y, z };
}

float xform_get_determinant(xform* m)
{
	return m->m[0][0] * m->m[1][1] * m->m[2][2]
		+ m->m[0][1] * m->m[1][2] * m->m[2][0]
		+ m->m[0][2] * m->m[1][0] * m->m[2][1]
		- m->m[2][0] * m->m[1][1] * m->m[0][2]
		- m->m[1][0] * m->m[0][1] * m->m[2][2]
		- m->m[0][0] * m->m[2][1] * m->m[1][2];
}