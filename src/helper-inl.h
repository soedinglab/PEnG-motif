/*
 * helper-inl.h
 *
 *  Created on: Jun 13, 2017
 *      Author: mmeier
 */

#ifndef SRC_HELPER_INL_H_
#define SRC_HELPER_INL_H_

template <typename E>
constexpr auto to_underlying(E e) noexcept
{
  return static_cast<std::underlying_type_t<E>>(e);
}


#endif /* SRC_HELPER_INL_H_ */
