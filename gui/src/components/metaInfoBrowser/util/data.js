export const updateList = (list, index, update) => ([
  ...list.slice(0, index),
  update(list[index]),
  ...list.slice(index + 1)
])

export const updateListState = (list, find, update, theDefault) => {
  let index = list.findIndex(find)
  if (index < 0) {
    return [...list, update(theDefault)]
  } else {
    return updateList(list, index, update)
  }
}
