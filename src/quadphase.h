#ifndef QUADTREEBD_H
#define QUADTREEBD_H

typedef struct vector2 {
  unsigned int x,y;
} vector2;

/* corner 0 is top left, going CCW */
typedef struct bbox2d {
  vector2 corners[4];
} bbox2d;

/*tl=top left br=bottom right*/
extern void bbox2d_set(bbox2d *bb, unsigned int tlx, unsigned int tly, unsigned int brx, unsigned int bry);
extern void vector2_set(vector2 *v, unsigned int x, unsigned int y);

typedef struct s_quadtree_node {
  bbox2d box;
  struct s_quadtree_node *childs[4];
  void *data;
} *quadtree_node;

#endif
