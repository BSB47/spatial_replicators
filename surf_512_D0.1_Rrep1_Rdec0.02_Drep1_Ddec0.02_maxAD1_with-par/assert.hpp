/* From Stroustrup special ed. p750 */

#ifndef ASSERTION
#define ASSERTION


template<class X,class A> inline void Assert(A assertion)
{
  if(!assertion) throw X();
}

class GeneralError {};

/* Look at p199 of Stroustrup special ed. (internal linking) */
#ifndef NDEBUG
#warning "Debugging mode is turned on"
namespace ASSERT{
  const bool ERROR_CHECK = true;
}
#else
namespace ASSERT{
  const bool ERROR_CHECK = false;
}
#endif

#endif
